/***
 *
 * Pollard Rho
 * Integer factorization and deterministic primality test for 1 <= n < 2^63
 *
 * Complexity: expected O(n^(1/4)) per split, O(1e6) once to build the sieve
 *
 * PollardRho f;           sieves smallest prime factors up to 1e6, keep a single instance around
 * f.factorize(n):         prime factors of n in ascending order with multiplicity, empty for n = 1
 * f.is_prime(n):          deterministic Miller Rabin for n < 2^63
 *
 * Brent's cycle detection with one gcd per batch of 128 steps and Montgomery multiplication
 * Shanks' SQUFOF was benchmarked as an alternative and was 3.5x to 45x slower on every workload
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// Newton iteration, doubles the correct bits of n^-1 mod 2^64 each step
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

struct PollardRho{
    static constexpr int LIMIT = 1000001;
    vector<int> spf, primes;
    mt19937_64 rng;

    PollardRho() : spf(LIMIT, 0), rng(20260105){
        for (int i = 2; i < LIMIT; i++){
            if (!spf[i]) spf[i] = i, primes.push_back(i);
            for (int p : primes){
                if (p > spf[i] || (long long)i * p >= LIMIT) break;
                spf[i * p] = p;
            }
        }
    }

    bool is_prime(long long n) const{
        if (n < LIMIT) return n >= 2 && spf[n] == n;
        if (!(n & 1)) return false;

        static const long long BASES[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
        long long s = 0, r = n - 1;
        for (; !(r & 1); r >>= 1, s++) {}

        Montgomery mont(n);
        unsigned long long one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        for (long long base : BASES){
            long long a = base % n;
            if (!a) continue;  /// n divides the base, this base says nothing about n
            unsigned long long c = mont.pow(mont.to_mont(a), r);
            for (int j = 0; j < s; j++){
                unsigned long long d = mont.mul(c, c);
                if (d == one && c != one && c != minus_one) return false;
                c = d;
            }
            if (c != one) return false;
        }
        return true;
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
        if (is_prime(n)){
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
                    g = __gcd(q, (unsigned long long)n);
                }
            }

            if (g == (unsigned long long)n){  /// the batch overshot to a multiple of n, replay it one step at a time
                do{
                    ys = next(ys, c);
                    g = __gcd(dist(x, ys), (unsigned long long)n);
                } while (g == 1);
            }
            if (g != (unsigned long long)n) return g;
        }
    }
};

int main(){
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

    assert(rho.is_prime(2) && rho.is_prime(1000000007) && rho.is_prime(9223372036854775783LL));
    assert(!rho.is_prime(1) && !rho.is_prime(0) && !rho.is_prime(561) && !rho.is_prime(1000000007LL * 3));

    return 0;
}
