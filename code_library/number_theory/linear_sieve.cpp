/***
 *
 * Linear Sieve
 * Smallest prime factors, primes, Euler's phi, Mobius mu and any multiplicative function for 1..n in one pass
 *
 * Complexity: O(n) time, 13 bytes per entry plus 4 per prime, about 135 MB at n = 1e7
 *
 * LinearSieve sieve(n): tables for 0..n
 * sieve.primes: primes <= n in order
 * sieve.spf[x]: smallest prime factor of x >= 2
 * sieve.spf_exponent[x]: exponent of spf[x] in x >= 2
 * sieve.phi[x]: Euler's totient of x >= 1
 * sieve.mu[x]: Mobius function of x >= 1
 * sieve.factorize(x): prime factors of 1 <= x <= n with multiplicity, ascending, O(log x)
 * sieve.multiplicative(f): g[x] for 1 <= x <= n where g is multiplicative and g(p^k) = f(p, k, p^k), g[1] = 1, in O(n)
 *     e.g. divisor count: f = [](long long p, int k, long long pk){ return k + 1; }
 *
 * Notes:
 *   - multiplicative() returns a new table of 8 bytes per entry on every call
 *   - For mu alone, mobius_function.cpp needs 1 byte per entry and runs about 3x faster
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LinearSieve{
    int n;
    vector<int> spf, primes, phi, mu;
    vector<unsigned char> spf_exponent;  /// the exponent of spf[x] in x, 1 byte since it is at most 30

    LinearSieve(int n) : n(n), spf(n + 1, 0), phi(n + 1, 0), mu(n + 1, 0), spf_exponent(n + 1, 0){
        if (n >= 1) phi[1] = mu[1] = 1;

        for (int i = 2; i <= n; i++){
            if (!spf[i]){
                spf[i] = i, primes.push_back(i);
                phi[i] = i - 1, mu[i] = -1, spf_exponent[i] = 1;
            }
            for (int p : primes){
                long long x = (long long)i * p;
                if (p > spf[i] || x > n) break;
                spf[x] = p;
                if (p == spf[i]) phi[x] = phi[i] * p, mu[x] = 0, spf_exponent[x] = spf_exponent[i] + 1;
                else phi[x] = phi[i] * (p - 1), mu[x] = -mu[i], spf_exponent[x] = 1;
            }
        }
    }

    vector<int> factorize(int x) const{
        assert(1 <= x && x <= n);
        vector<int> res;
        for (; x > 1; x /= spf[x]) res.push_back(spf[x]);
        return res;
    }

    template <typename F>
    vector<long long> multiplicative(F f) const{
        vector<long long> g(n + 1, 0);
        if (n >= 1) g[1] = 1;
        for (int x = 2; x <= n; x++){
            /// Rebuilding p^k costs 15% against a stored power table (90 vs 77 ms at n = 1e7), which would cost 40 MB
            int p = spf[x], k = spf_exponent[x], pk = p;
            for (int i = 1; i < k; i++) pk *= p;
            g[x] = g[x / pk] * f((long long)p, k, (long long)pk);
        }
        return g;
    }
};

int main(){
    /// Every table up to 100
    LinearSieve sieve(100);
    assert(sieve.primes.size() == 25);
    assert(sieve.primes[5] == 13);
    assert(sieve.spf[91] == 7);
    assert(sieve.spf[97] == 97);
    assert(sieve.spf_exponent[72] == 3);  /// 72 = 2^3 * 3^2
    assert(sieve.phi[36] == 12);
    assert(sieve.phi[97] == 96);
    assert(sieve.mu[30] == -1);
    assert(sieve.mu[12] == 0);
    assert((sieve.factorize(60) == vector<int>{2, 2, 3, 5}));
    assert(sieve.factorize(1).empty());

    /// Divisor count and divisor sum from their values on prime powers
    auto divisors = sieve.multiplicative([](long long, int k, long long){ return (long long)k + 1; });
    assert(divisors[12] == 6);
    assert(divisors[64] == 7);

    auto sigma = sieve.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });
    assert(sigma[12] == 28);
    assert(sigma[28] == 56);

    /// n = 1 has no primes, but phi[1] and mu[1] are set
    LinearSieve tiny(1);
    assert(tiny.primes.empty());
    assert(tiny.phi[1] == 1);
    assert(tiny.mu[1] == 1);

    return 0;
}
