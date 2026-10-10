/***
 *
 * Stirling Numbers
 * A whole row of Stirling numbers of the first or second kind modulo m, using exact NTT multiplication
 *
 * Complexity: O(n log n) for stirling_second and for stirling_first with a prime m > n,
 *             O(n log^2 n) for stirling_first with any other m
 *
 * stirling_first(n, m): c(n, k) for k = 0..n, the unsigned numbers counting permutations of n with k cycles
 *     coefficients of P_n(x) = x (x + 1) ... (x + n - 1), any modulus 1 <= m <= 2^30, n <= 2^22
 *     the signed numbers are s(n, k) = (-1)^(n - k) c(n, k)
 *     a prime m > n doubles with P_2k(x) = P_k(x) P_k(x + k), the Taylor shift dividing by factorials,
 *     any other m multiplies halves of the product instead, since factorials need not be invertible
 * stirling_second(n, p): S(n, k) for k = 0..n, the ways to split n labelled items into k non-empty groups
 *     S(n, k) = sum over i of (-1)^i / i! * (k - i)^n / (k - i)!, needs a prime p with n < p <= 2^30, n <= 2^22
 *
 * Products modulo 998244353 run one NTT, any other modulus runs three NTTs modulo 167772161, 469762049
 * and 754974721 joined by Garner's CRT, their product ~5.95e25 exceeds 2^23 * (2^30)^2, so every coefficient
 * is exact before the reduction modulo m
 * Holds no mutable global state, the roots are computed per product
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace stirling_ntt{
    const unsigned MOD = 998244353, P1 = 167772161, P2 = 469762049, P3 = 754974721;

    long long power(long long b, long long e, long long mod){
        long long res = 1 % mod;
        for (b %= mod; e; e >>= 1, b = b * b % mod){
            if (e & 1) res = res * b % mod;
        }
        return res;
    }

    /// Montgomery reduction with R = 2^32, x / R modulo P in [0, 2P) for any x < P R, replacing the slower x % P
    template<unsigned P>
    unsigned redc(unsigned long long x){
        constexpr unsigned NEG_INV = [](){
            unsigned y = P;
            for (int i = 0; i < 4; i++) y *= 2 - P * y;
            return -y;
        }();
        return (x + (unsigned long long)((unsigned)x * NEG_INV) * P) >> 32;
    }

    /// rt[k + j] = w^j R mod P for every power of two k < n, where w is a primitive 2k-th root of unity, or its inverse
    template<unsigned P, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), 1);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (P - 1) / (2 * k), P);
            if (invert) z = power(z, P - 2, P);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % P : rt[i >> 1];
        }

        for (auto& w : rt) w = ((unsigned long long)w << 32) % P;
        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, values stay in [0, 2P)
    /// P < 2^30 keeps u - v + 2P below 2^32 and every redc input below P R
    template<unsigned P>
    void forward(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();
        for (int k = n / 2; k >= 1; k >>= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = a[i + j + k];
                    a[i + j] = u + v >= 2 * P ? u + v - 2 * P : u + v;
                    a[i + j + k] = redc<P>((unsigned long long)(u - v + 2 * P) * rt[j + k]);
                }
            }
        }
    }

    /// Decimation in time with inverse roots, bit-reversed order in, natural order out scaled by n, values stay in [0, 2P)
    template<unsigned P>
    void inverse(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = redc<P>((unsigned long long)a[i + j + k] * rt[j + k]);
                    a[i + j] = u + v >= 2 * P ? u + v - 2 * P : u + v;
                    a[i + j + k] = u + 2 * P - v >= 2 * P ? u - v : u + 2 * P - v;
                }
            }
        }
    }

    /// Cyclic convolution of length n modulo P, inputs longer than n fold onto index i mod n
    template<unsigned P, unsigned ROOT>
    vector<unsigned> cyclic(const vector<long long>& a, const vector<long long>& b, int n){
        vector<unsigned> x(n, 0), y(n, 0);
        for (size_t i = 0; i < a.size(); i++) x[i & (n - 1)] = (x[i & (n - 1)] + a[i]) % P;
        for (size_t i = 0; i < b.size(); i++) y[i & (n - 1)] = (y[i & (n - 1)] + b[i]) % P;

        auto rt = roots<P, ROOT>(n, false);
        forward<P>(x, rt), forward<P>(y, rt);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * y[i] % P;

        rt = roots<P, ROOT>(n, true);
        inverse<P>(x, rt);

        unsigned long long inv = power(n, P - 2, P);
        for (auto& v : x) v = v * inv % P;
        return x;
    }

    /// Cyclic convolution of length n modulo m, for values in [0, m)
    vector<long long> cyclic_mod(const vector<long long>& a, const vector<long long>& b, int n, long long m){
        if (m == MOD){
            auto c = cyclic<MOD, 3>(a, b, n);
            return vector<long long>(c.begin(), c.end());
        }

        const unsigned long long INV1 = power(P1, P2 - 2, P2), INV12 = power((long long)P1 * P2 % P3, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (long long)P1 * P2 % m;
        auto c1 = cyclic<P1, 3>(a, b, n), c2 = cyclic<P2, 3>(a, b, n), c3 = cyclic<P3, 11>(a, b, n);

        /// x1 + x2 P1 + x3 P1 P2 is the true coefficient, reduced term by term below 2^62
        vector<long long> res(n);
        for (int i = 0; i < n; i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
        }
        return res;
    }

    /// a * b modulo m, for values in [0, m)
    vector<long long> multiply(const vector<long long>& a, const vector<long long>& b, long long m){
        if (a.empty() || b.empty()) return {};
        int need = a.size() + b.size() - 1, n = 1;
        while (n < need - 1) n <<= 1;

        /// A product one longer than n wraps only its top coefficient a.back() * b.back() onto index 0
        auto res = cyclic_mod(a, b, n, m);
        res.resize(need);
        if (need > n){
            long long top = a.back() * b.back() % m;
            res[0] = (res[0] - top + m) % m, res[n] = top;
        }
        return res;
    }

    bool is_prime(long long x){
        if (x < 2) return false;
        for (long long d = 2; d * d <= x; d++){
            if (x % d == 0) return false;
        }
        return true;
    }

    /// i! and 1 / i! modulo a prime p > n for i = 0..n
    void factorials(int n, long long p, vector<long long>& fact, vector<long long>& inv_fact){
        fact.assign(n + 1, 1), inv_fact.assign(n + 1, 1);
        for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % p;
        inv_fact[n] = power(fact[n], p - 2, p);
        for (int i = n; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % p;
    }
}

/// g(x + c) modulo a prime p, fact and inv_fact cover 0..deg g
vector<long long> taylor_shift(const vector<long long>& g, long long c, long long p, const vector<long long>& fact, const vector<long long>& inv_fact){
    int d = g.size() - 1;
    vector<long long> a(d + 1), b(d + 1);
    long long c_power = 1;
    for (int i = 0; i <= d; i++){
        a[i] = g[i] * fact[i] % p;
        b[d - i] = c_power * inv_fact[i] % p;
        c_power = c_power * c % p;
    }

    /// The j-th coefficient is (sum over i of a[i] c^(i - j) / (i - j)!) / j!, index d + j of a * b
    auto conv = stirling_ntt::multiply(a, b, p);
    vector<long long> res(d + 1);
    for (int j = 0; j <= d; j++) res[j] = conv[d + j] * inv_fact[j] % p;
    return res;
}

/// Product of (x + i) for l <= i < r modulo m
vector<long long> rising_product(int l, int r, long long m){
    if (r - l <= 32){
        vector<long long> res = {1 % m};
        for (int i = l; i < r; i++){
            res.push_back(0);
            for (int k = res.size() - 1; k > 0; k--) res[k] = (res[k - 1] + res[k] * (i % m)) % m;
            res[0] = res[0] * (i % m) % m;
        }
        return res;
    }

    int mid = (l + r) / 2;
    return stirling_ntt::multiply(rising_product(l, mid, m), rising_product(mid, r, m), m);
}

vector<long long> stirling_first(int n, long long m){
    assert(0 <= n && n <= (1 << 22) && 1 <= m && m <= (1LL << 30));
    if (n == 0 || n >= m || !stirling_ntt::is_prime(m)) return rising_product(0, n, m);

    vector<long long> fact, inv_fact;
    stirling_ntt::factorials(n / 2, m, fact, inv_fact);

    /// Walk the bits of n from the top, holding res = P_cur: doubling turns it into P_2cur, a set bit appends x + cur
    vector<long long> res = {1};
    int cur = 0;
    for (int bit = __lg(n); bit >= 0; bit--){
        if (cur){
            res = stirling_ntt::multiply(res, taylor_shift(res, cur, m, fact, inv_fact), m);
            cur *= 2;
        }

        if (n >> bit & 1){
            res.push_back(0);
            for (int k = cur + 1; k > 0; k--) res[k] = (res[k - 1] + res[k] * cur) % m;
            res[0] = res[0] * cur % m;
            cur++;
        }
    }
    return res;
}

vector<long long> stirling_second(int n, long long p){
    assert(0 <= n && n <= (1 << 22) && 2 <= p && n < p && p <= (1LL << 30) && stirling_ntt::is_prime(p));

    vector<long long> fact, inv_fact, x(n + 1), y(n + 1);
    stirling_ntt::factorials(n, p, fact, inv_fact);
    for (int i = 0; i <= n; i++){
        x[i] = stirling_ntt::power(i, n, p) * inv_fact[i] % p;
        y[i] = i % 2 ? (p - inv_fact[i]) % p : inv_fact[i];
    }

    vector<long long> res = stirling_ntt::multiply(x, y, p);
    res.resize(n + 1);
    return res;
}

int main(){
    assert((stirling_first(0, 1000000007) == vector<long long>{1}));
    assert((stirling_first(1, 1000000007) == vector<long long>{0, 1}));
    assert((stirling_first(4, 1000000007) == vector<long long>{0, 6, 11, 6, 1}));
    assert((stirling_first(5, 1000000007) == vector<long long>{0, 24, 50, 35, 10, 1}));
    assert((stirling_first(5, 998244353) == vector<long long>{0, 24, 50, 35, 10, 1}));
    assert((stirling_first(5, 7) == vector<long long>{0, 3, 1, 0, 3, 1}));
    assert((stirling_first(7, 7) == vector<long long>{0, 6, 0, 0, 0, 0, 0, 1}));
    assert((stirling_first(5, 12) == vector<long long>{0, 0, 2, 11, 10, 1}));
    assert((stirling_first(3, 1) == vector<long long>{0, 0, 0, 0}));

    assert((stirling_second(0, 1000000007) == vector<long long>{1}));
    assert((stirling_second(4, 1000000007) == vector<long long>{0, 1, 7, 6, 1}));
    assert((stirling_second(5, 1000000007) == vector<long long>{0, 1, 15, 25, 10, 1}));
    assert((stirling_second(5, 998244353) == vector<long long>{0, 1, 15, 25, 10, 1}));
    assert((stirling_second(5, 7) == vector<long long>{0, 1, 1, 4, 3, 1}));

    return 0;
}
