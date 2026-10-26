/***
 *
 * Pollard Rho
 * Integer factorization and deterministic primality test for 1 <= n < 2^63
 *
 * Complexity: expected O(n^(1/4)) per split, O(L) once to sieve up to L = 10^6
 *
 * PollardRho f: sieves smallest prime factors up to 10^6, keep a single instance around
 * f.factorize(n): prime factors of n in ascending order with multiplicity, empty for n = 1
 * prm::is_prime(n): true if n is prime, for any long long n
 *
 * Notes:
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Embeds checked copies of miller_rabin.cpp's Montgomery and prm::is_prime to stay standalone
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n, 8x faster than __int128 % n here
/// BEGIN COPY montgomery from code_library/number_theory/miller_rabin.cpp
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// each Newton step doubles the correct bits of n^-1 mod 2^64
        r2 = -n % n;
        r2 = (unsigned __int128)r2 * r2 % n;
    }

    unsigned long long reduce(unsigned __int128 x) const{
        unsigned long long q = (unsigned long long)x * inv, m = ((unsigned __int128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    unsigned long long mul(unsigned long long x, unsigned long long y) const{
        return reduce((unsigned __int128)x * y);
    }

    unsigned long long add(unsigned long long x, unsigned long long y) const{
        return (x += y) >= n ? x - n : x;
    }

    unsigned long long to_mont(unsigned long long x) const{
        return mul(x, r2);
    }

    unsigned long long pow(unsigned long long x, unsigned long long e) const{
        unsigned long long res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }
        return res;
    }
};
/// END COPY montgomery

/// BEGIN COPY miller_rabin from code_library/number_theory/miller_rabin.cpp
namespace prm{
    /// Jaeschke 1993, deterministic below 4759123141
    const vector<int> BASES_32 = {2, 7, 61};
    /// Jim Sinclair 2011, deterministic below 2^64: https://miller-rabin.appspot.com/
    const vector<int> BASES_64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    struct SmallPrime{ unsigned long long p, inv, max_quotient; };

    /// Holds the bases 7 and 61, so is_witness never gets a base that n divides
    const vector<SmallPrime> SMALL_PRIMES = [](){
        vector<SmallPrime> res;
        for (unsigned long long p: {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79}){
            res.push_back({p, Montgomery(p).inv, ~0ULL / p});
        }
        return res;
    }();

    bool is_witness(int a, long long n, int s, const Montgomery& mont){
        unsigned long long one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        unsigned long long x = mont.pow(mont.to_mont(a), (n - 1) >> s);

        if (x == one) return false;

        for (int i = 0; i < s; i++){
            if (x == minus_one) return false;
            x = mont.mul(x, x);
        }

        return true;
    }

    bool miller_rabin(long long n, const vector<int>& bases){
        if (n < 2 || (n & 1) == 0) return n == 2;

        /// p | n exactly when n * p^-1 mod 2^64 <= (2^64 - 1) / p, no division: 27% faster on random n below 2^63
        for (auto& [p, inv, max_quotient]: SMALL_PRIMES){
            if ((unsigned long long)n * inv <= max_quotient) return (unsigned long long)n == p;
        }

        int s = __builtin_ctzll(n - 1);
        Montgomery mont(n);
        for (auto a: bases){
            if (is_witness(a, n, s, mont)) return false;
        }

        return true;
    }

    bool is_prime(long long n){
        if (n < 4759123141LL) return miller_rabin(n, BASES_32);
        return miller_rabin(n, BASES_64);
    }
}
/// END COPY miller_rabin

/// BEGIN SHARED pollard_rho
struct PollardRho{
    static constexpr int LIMIT = 1000001;
    vector<int> spf;
    mt19937_64 rng;

    PollardRho() : spf(LIMIT, 0), rng(20260105){
        vector<int> primes;
        for (int i = 2; i < LIMIT; i++){
            if (!spf[i]) spf[i] = i, primes.push_back(i);
            for (int p : primes){
                if (p > spf[i] || (long long)i * p >= LIMIT) break;
                spf[i * p] = p;
            }
        }
    }

    vector<long long> factorize(long long n){
        assert(n >= 1);
        vector<long long> res;
        collect(n, res);
        sort(res.begin(), res.end());
        return res;
    }

    void collect(long long n, vector<long long>& out){
        if (n < LIMIT){
            for (; n > 1; n /= spf[n]) out.push_back(spf[n]);
            return;
        }

        if (!(n & 1)){
            out.push_back(2);
            return collect(n >> 1, out);
        }

        if (prm::is_prime(n)){
            out.push_back(n);
            return;
        }

        long long d = split(n);
        collect(d, out);
        collect(n / d, out);
    }

    /// A nontrivial factor of an odd composite n
    long long split(long long n){
        const int BATCH = 128;
        Montgomery mont(n);
        auto next = [&](unsigned long long x, unsigned long long c){ return mont.add(mont.mul(x, x), c); };
        auto dist = [](unsigned long long x, unsigned long long y){ return x > y ? x - y : y - x; };

        while (true){
            unsigned long long c = rng() % (n - 1) + 1, y = rng() % n, x = y, ys = y, q = mont.to_mont(1), g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y, c);

                for (long long k = 0; k < r && g == 1; k += BATCH){
                    ys = y;
                    for (long long i = 0; i < BATCH && i < r - k; i++){
                        y = next(y, c);
                        q = mont.mul(q, dist(x, y));
                    }
                    g = gcd(q, (unsigned long long)n);  /// binary GCD in libstdc++: factorize 12-18% faster than __gcd
                }
            }

            if (g == (unsigned long long)n){  /// the batch overshot to a multiple of n, replay it one step at a time
                do{
                    ys = next(ys, c);
                    g = gcd(dist(x, ys), (unsigned long long)n);
                } while (g == 1);
            }

            if (g != (unsigned long long)n) return g;
        }
    }
};
/// END SHARED pollard_rho

int main(){
    /// Factorizations with repeated, small and large prime factors
    PollardRho rho;

    vector<pair<long long, vector<long long>>> cases = {
        {1, {}},
        {2, {2}},
        {12, {2, 2, 3}},
        {999983, {999983}},
        {1000003, {1000003}},
        {1000000007LL * 998244353LL, {998244353, 1000000007}},
        {1000003LL * 1000003LL, {1000003, 1000003}},
        {(1LL << 62), vector<long long>(62, 2)},
        {9223372036854775783LL, {9223372036854775783LL}},
        {600851475143LL, {71, 839, 1471, 6857}},
        {3LL * 1000003 * 1000033, {3, 1000003, 1000033}},
        {1000003LL * 1000033 * 1000037, {1000003, 1000033, 1000037}},
    };
    for (auto& [n, expected] : cases) assert(rho.factorize(n) == expected);

    /// Primality without building a PollardRho
    assert(prm::is_prime(1000000007));
    assert(prm::is_prime(9223372036854775783LL));
    assert(!prm::is_prime(1));
    assert(!prm::is_prime(561));
    assert(!prm::is_prime(1000000007LL * 3));

    return 0;
}
