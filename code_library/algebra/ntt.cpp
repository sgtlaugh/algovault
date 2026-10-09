/***
 *
 * Number Theoretic Transform
 * Polynomial multiplication modulo an NTT-friendly prime, modulo any m < 2^31, or exact in 64 bits
 *
 * Complexity: O(n log n) for a product of total length n
 *
 * ntt::multiply<MOD, ROOT>(a, b): a * b modulo MOD, default MOD = 998244353 and ROOT = 3
 *     MOD must be a prime below 2^31 and ROOT a quadratic non-residue modulo MOD (any primitive root is one)
 *     the result length |a| + |b| - 1, rounded up to a power of two, must divide MOD - 1:
 *     at most 2^23 for 998244353, 2^20 for 7340033 = 7 * 2^20 + 1
 * ntt::mod_multiply(a, b, m): a * b modulo any 1 <= m < 2^31, prime or not, result length at most 2^24
 * ntt::exact_multiply(a, b): the exact a * b, correct whenever every true coefficient lies in [-2^63, 2^63)
 *     result length at most 2^24, inputs may be any long long
 * ntt::crt_multiply(a, b): the core of the two above, every true coefficient of a * b modulo M as an unsigned __int128 in [0, M)
 *     result length at most 2^24
 *
 * crt_multiply runs three NTTs modulo 167772161, 469762049 and 754974721 and combines them with Garner's CRT,
 * their product M ~ 5.95e25 exceeds 2^23 * (2^31 - 1)^2, so every mod_multiply coefficient is recovered exactly
 * Inputs may be negative or >= the modulus, they are reduced first, results lie in [0, modulus)
 * An empty input gives an empty result
 * Holds no mutable global state, the roots are computed per call in O(n)
 * Requires __int128 (64-bit GCC or Clang) for mod_multiply and exact_multiply
 *
 * Example:
 *     ntt::multiply({1, 2, 3}, {4, 5}) == {4, 13, 22, 15}
 *     ntt::mod_multiply({1000000006}, {1000000006}, 1000000007) == {1}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace ntt{
    const unsigned P1 = 167772161, P2 = 469762049, P3 = 754974721;

    constexpr unsigned long long power(unsigned long long b, unsigned long long e, unsigned long long mod){
        unsigned long long res = 1;
        for (b %= mod; e; e >>= 1, b = b * b % mod){
            if (e & 1) res = res * b % mod;
        }
        return res;
    }

    constexpr bool is_prime(unsigned n){
        if (n < 2) return false;
        for (unsigned long long d = 2; d * d <= n; d++){
            if (n % d == 0) return false;
        }
        return true;
    }

    /// rt[k + j] = w^j for every power of two k < n, where w is a primitive 2k-th root of unity
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n){
        vector<unsigned> rt(max(n, 2), 1);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % MOD : rt[i >> 1];
        }
        return rt;
    }

    /// Forward transform in place, a.size() must be a power of two no larger than rt.size()
    template<unsigned MOD>
    void transform(vector<unsigned>& a, const vector<unsigned>& rt){
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++){
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }

        /// MOD < 2^31 keeps u + v below 2^32
        for (int k = 1; k < n; k <<= 1){
            for (int i = 0; i < n; i += 2 * k){
                for (int j = 0; j < k; j++){
                    unsigned u = a[i + j], v = (unsigned long long)rt[j + k] * a[i + j + k] % MOD;
                    a[i + j] = u + v >= MOD ? u + v - MOD : u + v;
                    a[i + j + k] = u >= v ? u - v : u + MOD - v;
                }
            }
        }
    }

    template<unsigned MOD = 998244353, unsigned ROOT = 3>
    vector<long long> multiply(const vector<long long>& a, const vector<long long>& b){
        static_assert(MOD < (1u << 31) && is_prime(MOD), "MOD must be a prime below 2^31");
        static_assert(power(ROOT, (MOD - 1) / 2, MOD) == MOD - 1, "ROOT must be a quadratic non-residue modulo MOD");
        if (a.empty() || b.empty()) return {};

        int s = a.size() + b.size() - 1, n = 1;
        while (n < s) n <<= 1;
        assert((MOD - 1) % n == 0);

        vector<unsigned> x(n, 0), y(n, 0);
        for (size_t i = 0; i < a.size(); i++) x[i] = (a[i] % (long long)MOD + MOD) % MOD;
        for (size_t i = 0; i < b.size(); i++) y[i] = (b[i] % (long long)MOD + MOD) % MOD;

        auto rt = roots<MOD, ROOT>(n);
        transform<MOD>(x, rt), transform<MOD>(y, rt);

        /// The forward transform with reversed indices 1..n-1 is the inverse transform scaled by n
        unsigned long long inv = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * y[i] % MOD * inv % MOD;
        reverse(x.begin() + 1, x.end());
        transform<MOD>(x, rt);

        return vector<long long>(x.begin(), x.begin() + s);
    }

    /// The true coefficients modulo M = P1 * P2 * P3, in [0, M)
    vector<unsigned __int128> crt_multiply(const vector<long long>& a, const vector<long long>& b){
        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);

        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        vector<unsigned __int128> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = x1 + (unsigned __int128)x2 * P1 + (unsigned __int128)x3 * P1 * P2;
        }

        return res;
    }

    vector<long long> mod_multiply(vector<long long> a, vector<long long> b, long long m){
        assert(1 <= m && m < (1LL << 31));
        for (auto& x : a) x = (x % m + m) % m;
        for (auto& x : b) x = (x % m + m) % m;

        auto c = crt_multiply(a, b);
        vector<long long> res(c.size());
        for (size_t i = 0; i < c.size(); i++) res[i] = c[i] % m;
        return res;
    }

    vector<long long> exact_multiply(const vector<long long>& a, const vector<long long>& b){
        const unsigned __int128 M = (unsigned __int128)P1 * P2 * P3;

        auto c = crt_multiply(a, b);
        vector<long long> res(c.size());
        for (size_t i = 0; i < c.size(); i++) res[i] = c[i] > M / 2 ? (long long)((__int128)c[i] - (__int128)M) : (long long)c[i];
        return res;
    }
}

int main(){
    using namespace ntt;

    const long long MOD = 998244353;

    assert((multiply({1, 2, 3}, {4, 5}) == vector<long long>{4, 13, 22, 15}));
    assert((multiply({5, 1, 2, 6, 9, 8}, {3, 9, 0, 2}) == vector<long long>{15, 48, 15, 46, 83, 109, 84, 18, 16}));
    assert((multiply({7}, {6}) == vector<long long>{42}));
    assert((multiply({}, {1, 2}) == vector<long long>{}));
    assert((multiply({MOD - 1, MOD - 1}, {MOD - 1}) == vector<long long>{1, 1}));
    assert((multiply({-1}, {2}) == vector<long long>{MOD - 2}));
    assert((multiply({MOD + 3}, {2 * MOD + 5}) == vector<long long>{15}));
    assert((multiply<17, 3>({16, 16}, {16, 16}) == vector<long long>{1, 2, 1}));
    assert((multiply<7340033, 3>({1, 1}, {1, 1, 1}) == vector<long long>{1, 2, 2, 1}));

    assert((mod_multiply({5, 1, 2, 6, 9, 8}, {3, 9, 0, 2}, 14) == vector<long long>{1, 6, 1, 4, 13, 11, 0, 4, 2}));
    assert((mod_multiply({1000000006, 1000000006}, {1000000006, 1000000006}, 1000000007) == vector<long long>{1, 2, 1}));
    assert((mod_multiply({2147483646}, {2147483646}, 2147483647) == vector<long long>{1}));
    assert((mod_multiply({-3, 4}, {5}, 7) == vector<long long>{6, 6}));
    assert((mod_multiply({123, 456}, {789}, 1) == vector<long long>{0, 0}));

    assert((exact_multiply({1000000000, 1000000000}, {1000000000, 1000000000}) == vector<long long>{1000000000000000000LL, 2000000000000000000LL, 1000000000000000000LL}));
    assert((exact_multiply({3037000499LL}, {3037000499LL}) == vector<long long>{9223372030926249001LL}));
    assert((exact_multiply({LLONG_MAX, LLONG_MIN}, {1}) == vector<long long>{LLONG_MAX, LLONG_MIN}));
    assert((exact_multiply({-2147483648LL}, {4294967296LL}) == vector<long long>{LLONG_MIN}));
    assert((exact_multiply({-3, 2}, {4, -1}) == vector<long long>{-12, 11, -2}));

    return 0;
}
