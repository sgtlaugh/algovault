/***
 *
 * Number Theoretic Transform
 * Polynomial multiplication modulo an NTT-friendly prime, modulo any m < 2^31, or exact in 64 bits
 *
 * Complexity: O(n log n) for a product of total length n
 *
 * ntt::multiply<MOD, ROOT>(a, b): a * b modulo MOD, default MOD = 998244353, ROOT found from MOD, the result length rounded up to a power of two must divide MOD - 1
 * ntt::mod_multiply(a, b, m): a * b modulo any 1 <= m < 2^31, prime or not, result length at most 2^24
 * ntt::exact_multiply(a, b): the exact a * b of any long long inputs, correct whenever every true coefficient lies in [-2^63, 2^63), result length at most 2^24
 * ntt::crt_multiply(a, b): every true coefficient of a * b modulo M = 167772161 * 469762049 * 754974721 ~ 5.95e25 as an unsigned __int128 in [0, M), result length at most 2^24
 *
 * The result length is |a| + |b| - 1, and an empty input gives an empty result
 * MOD must be a prime below 2^31; an explicit ROOT must be a quadratic non-residue modulo MOD, which any primitive root is
 * multiply allows result lengths up to 2^23 modulo 998244353 and 2^20 modulo 7340033 = 7 * 2^20 + 1
 * Inputs may be negative or >= the modulus, multiply and mod_multiply return values in [0, modulus)
 * Requires __int128 (64-bit GCC or Clang) for exact_multiply and crt_multiply
 * Use fft.cpp for exact products of small values like 0..100, about 3x faster there; any other mod runs at a similar speed in both
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

namespace ntt{
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

    /// Values live in [0, BOUND): 2 MOD below 2^30 saves the final reduction, MOD above, so u - v + BOUND stays below 2^32
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

    /// rt[k + j] = w^j R modulo MOD for every power of two k < n, where w is a primitive 2k-th root of unity, or its inverse
    template<unsigned MOD, unsigned ROOT>
    vector<unsigned> roots(int n, bool invert){
        vector<unsigned> rt(max(n, 2), 1);
        for (int k = 2; k < n; k <<= 1){
            unsigned long long z = power(ROOT, (MOD - 1) / (2 * k), MOD);
            if (invert) z = power(z, MOD - 2, MOD);
            for (int i = k; i < 2 * k; i++) rt[i] = i & 1 ? rt[i >> 1] * z % MOD : rt[i >> 1];
        }

        for (auto& w : rt) w = ((unsigned long long)w << 32) % MOD;
        return rt;
    }

    /// Decimation in frequency, natural order in, bit-reversed order out, a.size() a power of two no larger than rt.size()
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

    /// Any quadratic non-residue g works as ROOT: g^((MOD - 1) / n) then has order exactly n for every power of two n dividing MOD - 1
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

        vector<unsigned> x(n, 0), y(n, 0);
        for (size_t i = 0; i < a.size(); i++) x[i] = (a[i] % (long long)MOD + MOD) % MOD;
        for (size_t i = 0; i < b.size(); i++) y[i] = (b[i] % (long long)MOD + MOD) % MOD;

        auto rt = roots<MOD, ROOT>(n, false);
        forward_transform<MOD>(x, rt);
        forward_transform<MOD>(y, rt);

        unsigned long long inv = power(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++) x[i] = (unsigned long long)x[i] * y[i] % MOD * inv % MOD;
        inverse_transform<MOD>(x, roots<MOD, ROOT>(n, true));

        vector<long long> res(s);
        for (int i = 0; i < s; i++) res[i] = x[i] >= MOD ? x[i] - MOD : x[i];
        return res;
    }

    vector<unsigned __int128> crt_multiply(const vector<long long>& a, const vector<long long>& b){
        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);

        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        vector<unsigned __int128> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            /// Garner: res = x1 + x2 * P1 + x3 * P1 * P2, each digit chosen so res matches c1, c2, c3 modulo P1, P2, P3
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

        constexpr unsigned long long INV1 = power(P1, P2 - 2, P2);
        constexpr unsigned long long INV12 = power((unsigned long long)P1 * P2, P3 - 2, P3);
        const unsigned long long M1 = P1 % m, M12 = (unsigned long long)P1 * P2 % m;
        auto c1 = multiply<P1, 3>(a, b), c2 = multiply<P2, 3>(a, b), c3 = multiply<P3, 11>(a, b);

        /// Every true coefficient is below 2^23 m^2 < M ~ 5.95e25, so x1 + x2 P1 + x3 P1 P2 is exact, reduced term by term below 2^62
        vector<long long> res(c1.size());
        for (size_t i = 0; i < res.size(); i++){
            unsigned long long x1 = c1[i];
            unsigned long long x2 = (c2[i] + P2 - x1) % P2 * INV1 % P2;
            unsigned long long x3 = (c3[i] + 2ULL * P3 - x1 - x2 * P1 % P3) % P3 * INV12 % P3;
            res[i] = (x1 + x2 * M1 + x3 * M12) % m;
        }
        return res;
    }

    vector<long long> exact_multiply(const vector<long long>& a, const vector<long long>& b){
        const unsigned __int128 M = (unsigned __int128)P1 * P2 * P3;

        auto c = crt_multiply(a, b);
        vector<long long> res(c.size());

        /// A negative true value v comes back as M + v, which lies above M / 2
        for (size_t i = 0; i < c.size(); i++) res[i] = c[i] > M / 2 ? (long long)((__int128)c[i] - (__int128)M) : (long long)c[i];
        return res;
    }
}

int main(){
    using namespace ntt;

    /// (1 + 2x + 3x^2)(4 + 5x) modulo 998244353, where -1 * 2 comes back as 998244353 - 2
    assert((multiply({1, 2, 3}, {4, 5}) == vector<long long>{4, 13, 22, 15}));
    assert((multiply({-1}, {2}) == vector<long long>{998244351}));

    /// Another NTT-friendly prime, 786433 = 3 * 2^18 + 1, where 3 is a quadratic residue so its root is found as 5
    assert((multiply<786433>({1, 1}, {1, 1, 1}) == vector<long long>{1, 2, 2, 1}));

    /// 1000000006 is -1 modulo 10^9 + 7, so this is (-1 - x)^2, and 14 is a composite modulus
    assert((mod_multiply({1000000006, 1000000006}, {1000000006, 1000000006}, 1000000007) == vector<long long>{1, 2, 1}));
    assert((mod_multiply({5, 1, 2, 6, 9, 8}, {3, 9, 0, 2}, 14) == vector<long long>{1, 6, 1, 4, 13, 11, 0, 4, 2}));

    /// Signed results, and (10^9 + 10^9 x)^2 = 10^18 + 2 * 10^18 x + 10^18 x^2, past any single prime
    assert((exact_multiply({-3, 2}, {4, -1}) == vector<long long>{-12, 11, -2}));
    assert((exact_multiply({1000000000, 1000000000}, {1000000000, 1000000000}) == vector<long long>{1000000000000000000LL, 2000000000000000000LL, 1000000000000000000LL}));

    /// 2^40 * 2^40 = 2^80 does not fit in 64 bits but lies below M
    assert((crt_multiply({1LL << 40}, {1LL << 40}) == vector<unsigned __int128>{(unsigned __int128)1 << 80}));

    return 0;
}
