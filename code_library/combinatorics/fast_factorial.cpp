/***
 *
 * Fast Factorial
 * n! modulo a prime p below 2^31 for any n, without a loop up to n
 *
 * Complexity, for the block size b = fact.block:
 *   - constructor: O((p / b) log b), which is O(sqrt(p) log p) for one query
 *   - fact(n): O(b)
 *   - memory: 12 p / b bytes while building, 4 p / b bytes kept
 *
 * FastFactorial fact(p, queries = 1): tables for a prime 2 <= p < 2^31, sized for about queries >= 1 calls of fact
 * fact(n): n! modulo p for any n >= 0, 0 when n >= p
 * fact.block: the largest power of two with b^2 + b / 2 < p, halved towards sqrt(128 p / queries) but not below 256
 *
 * Notes:
 *   - p must be prime, which is not checked; any prime works, 10^9 + 7 as well as 998244353
 *   - queries only trades build time against call time, any number of calls stays correct
 *   - For n below about 10^7, the O(n) tables of combinatorics.cpp are simpler and answer in O(1)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace fast_factorial_ntt{
    /// BEGIN COPY ntt_multiply from code_library/algebra/ntt.cpp
    const unsigned P1 = 167772161, P2 = 469762049, P3 = 754974721;

    constexpr bool is_prime(unsigned n){
        if (n < 2) return false;
        for (unsigned long long d = 2; d * d <= n; d++){
            if (n % d == 0) return false;
        }
        return true;
    }

    /// BEGIN SHARED ntt_core
    constexpr unsigned long long power(unsigned long long b, unsigned long long e, unsigned long long mod){
        unsigned long long res = 1;
        b %= mod;
        while (e){
            if (e & 1) res = res * b % mod;
            b = b * b % mod;
            e >>= 1;
        }
        return res;
    }

    /// Newton's step y (2 - p y) doubles the correct low bits of 1 / p modulo 2^32, y = p starts with 3 of them
    constexpr unsigned inverse_mod_2_32(unsigned p, unsigned y, int steps){
        return steps ? inverse_mod_2_32(p, y * (2 - p * y), steps - 1) : y;
    }

    /// Values live in [0, BOUND), 2 MOD below 2^30 to skip the final reduction, else MOD, so u - v + BOUND < 2^32
    template<unsigned MOD>
    constexpr unsigned bound(){
        return MOD < (1u << 30) ? 2 * MOD : MOD;
    }

    /// Montgomery reduction, x / 2^32 modulo MOD for x < MOD 2^32, in [0, bound<MOD>()), replacing the slower x % MOD
    template<unsigned MOD>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = 0u - inverse_mod_2_32(MOD, MOD, 4);
        unsigned r = (x + (unsigned long long)((unsigned)x * NEG_INV) * MOD) >> 32;
        return MOD < (1u << 30) || r < MOD ? r : r - MOD;
    }

    /// rt[k + j] = w^j R modulo MOD for each power of two k < n, w a primitive 2k-th root of unity or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), ((unsigned long long)1 << 32) % MOD);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            z = (z << 32) % MOD;
            for (int i = k; i < 2 * k; i++){
                unsigned w = redc<MOD>(rt[i >> 1] * z);
                rt[i] = i & 1 ? (w >= MOD ? w - MOD : w) : rt[i >> 1];
            }
        }
        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two at most rt.size()
    template<unsigned MOD>
    void forward_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = redc<MOD>((unsigned long long)(u - v + BOUND) * rt[j + k]);
                }
            }
        }
    }

    /// Decimation in time with the inverse roots, bit-reversed order in, natural order out, scaled by n
    template<unsigned MOD>
    void inverse_transform(vector<unsigned>& a, const vector<unsigned>& rt){
        const unsigned BOUND = bound<MOD>();
        int n = a.size();
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<MOD>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= BOUND ? u + v - BOUND : u + v;
                    a[i + j + k] = u >= v ? u - v : u + BOUND - v;
                }
            }
        }
    }
    /// END SHARED ntt_core

    /// Any quadratic non-residue g works as ROOT: g^((MOD - 1) / n) has order n for each power of two n | MOD - 1
    constexpr unsigned non_residue(unsigned mod){
        unsigned g = 2;
        while (power(g, (mod - 1) / 2, mod) != mod - 1) g++;
        return g;
    }

    template<unsigned MOD = 998244353, unsigned ROOT = non_residue(MOD)>
    vector<long long> multiply(const vector<long long>& a, const vector<long long>& b){
        static_assert(MOD < (1u << 31) && is_prime(MOD), "MOD must be a prime below 2^31");
        static_assert(power(ROOT, (MOD - 1) / 2, MOD) == MOD - 1, "ROOT must be a quadratic non-residue modulo MOD");
        if (a.empty() || b.empty()) return {};

        int s = a.size() + b.size() - 1, n = 1;
        while (n < s) n <<= 1;
        assert((MOD - 1) % n == 0);

        /// Each level's roots sit at fixed positions, so one table grown to the largest n serves every smaller call
        static vector<unsigned> rt, irt;
        if ((int)rt.size() < n) rt = roots<MOD, ROOT>(n, false), irt = roots<MOD, ROOT>(n, true);

        bool same = a == b;
        vector<unsigned> x(n, 0), y;
        for (size_t i = 0; i < a.size(); i++) x[i] = (a[i] % (long long)MOD + MOD) % MOD;
        forward_transform<MOD>(x, rt);
        if (!same){
            y.assign(n, 0);
            for (size_t i = 0; i < b.size(); i++) y[i] = (b[i] % (long long)MOD + MOD) % MOD;
            forward_transform<MOD>(y, rt);
        }

        const vector<unsigned>& z = same ? x : y;
        unsigned long long inv = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * z[i] % MOD * inv % MOD;
        inverse_transform<MOD>(x, irt);

        vector<long long> res(s);
        for (int i = 0; i < s; i++) res[i] = x[i] >= MOD ? x[i] - MOD : x[i];
        return res;
    }
    /// END COPY ntt_multiply

    /// BEGIN COPY ntt_mod_multiply from code_library/algebra/ntt.cpp
    vector<long long> mod_multiply(vector<long long> a, vector<long long> b, long long m){
        assert(1 <= m && m < (1LL << 31));
        for (auto& x : a) x = (x % m + m) % m;
        for (auto& x : b) x = (x % m + m) % m;

        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (unsigned long long)P1 * P2 % m;
        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        /// True coefficients < 2^23 m^2 < M ~ 5.95e25 make x1 + x2 P1 + x3 P1 P2 exact, reduced below 2^62 per term
        vector<long long> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
        }

        return res;
    }
    /// END COPY ntt_mod_multiply
}

struct FastFactorial{
    static const int MIN_BLOCK = 256, COST = 128;

    long long p;
    int block;
    /// table[q] = (q block)! modulo p
    vector<int> table;

    FastFactorial(long long p, long long queries = 1) : p(p){
        assert(2 <= p && p < (1LL << 31) && queries >= 1);

        /// Shifting by d / block lands a sample point on 0 modulo p unless block^2 + block / 2 < p
        block = 1;
        while (4LL * block * block + block < p) block *= 2;
        /// A table entry costs about COST multiplications, a query up to block, so many queries shrink the block
        while (block > MIN_BLOCK && (p / block) * COST / queries < block) block /= 2;

        int count = (p - 1) / block + 1;
        vector<long long> g = block_products(count - 1);
        table.assign(count, 1);
        for (int q = 1; q < count; q++) table[q] = table[q - 1] * g[q - 1] % p;
    }

    long long operator()(long long n){
        assert(n >= 0);
        if (n >= p) return 0;

        long long q = n / block, res = table[q];
        for (long long i = q * block + 1; i <= n; i++) res = res * i % p;
        return res;
    }

    /// g[i] = (i block + 1) (i block + 2) ... (i block + block) modulo p for i = 0..k - 1
    vector<long long> block_products(int k){
        long long inv_block = fast_factorial_ntt::power(block, p - 2, p);

        /// h[x] = g_d(x) = (block x + 1) ... (block x + d) for x = 0..d, and g_2d(x) = g_d(x) g_d(x + d / block)
        vector<long long> h = {1, (block + 1) % p};
        for (int d = 1; d < block; d *= 2){
            vector<long long> right = shift(h, d + 1, d - 1), moved = shift(h, d * inv_block % p, 2 * d);
            h.insert(h.end(), right.begin(), right.end());
            for (int x = 0; x <= 2 * d; x++) h[x] = h[x] * moved[x] % p;
        }

        vector<long long> base = h;
        while ((int)h.size() < k){
            vector<long long> more = shift(base, h.size(), block);
            h.insert(h.end(), more.begin(), more.end());
        }
        h.resize(k);
        return h;
    }

    /// f(a), ..., f(a + m) for the degree d polynomial f with f(i) = h[i], needs a - d, ..., a + m nonzero modulo p
    vector<long long> shift(const vector<long long>& h, long long a, int m){
        int d = h.size() - 1, len = d + m + 1;

        long long fact = 1;
        for (int i = 1; i <= d; i++) fact = fact * i % p;
        vector<long long> inv_fact(d + 1);
        inv_fact[d] = fast_factorial_ntt::power(fact, p - 2, p);
        for (int i = d; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % p;

        /// Lagrange: f(a + k) = prod_j (a + k - j) * sum_i c[i] / (a + k - i), the sum a convolution of c with 1 / x
        vector<long long> c(d + 1);
        for (int i = 0; i <= d; i++){
            c[i] = h[i] * inv_fact[i] % p * inv_fact[d - i] % p;
            if ((d - i) & 1) c[i] = (p - c[i]) % p;
        }

        vector<long long> x(len), pre(len + 1, 1), inv_pre(len + 1), inv_x(len);
        for (int t = 0; t < len; t++){
            x[t] = ((a - d + t) % p + p) % p;
            pre[t + 1] = pre[t] * x[t] % p;
        }
        inv_pre[len] = fast_factorial_ntt::power(pre[len], p - 2, p);
        for (int t = len - 1; t >= 0; t--){
            inv_pre[t] = inv_pre[t + 1] * x[t] % p;
            inv_x[t] = inv_pre[t + 1] * pre[t] % p;
        }

        vector<long long> conv = fast_factorial_ntt::mod_multiply(c, inv_x, p), res(m + 1);
        for (int k = 0; k <= m; k++) res[k] = conv[d + k] * pre[k + d + 1] % p * inv_pre[k] % p;
        return res;
    }
};

int main(){
    /// 13! = 6227020800 = 6 (10^9 + 7) + 227020758, and Wilson's theorem gives (p - 1)! = -1 and (p - 2)! = 1
    FastFactorial fact(1000000007);
    assert(fact(0) == 1);
    assert(fact(13) == 227020758);
    assert(fact(1000000006) == 1000000006);
    assert(fact(1000000005) == 1);
    assert(fact(1000000007) == 0);

    /// Sized for many calls modulo 998244353, which shrinks the block and so the cost of each call
    FastFactorial many(998244353, 1000);
    assert(many(10) == 3628800);

    return 0;
}
