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
 *   - n up to about 1e8 by memory, 1.3 GB; n = 2^31 - 1 overflows n + 1
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

#define MAX 500000

int main(){
    int t, n;

    LinearSieve sieve(MAX);
    vector<long long> sigma = sieve.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });

    scanf("%d", &t);
    while (t--){
        scanf("%d", &n);
        printf("%lld\n", sigma[n] - n);
    }

    return 0;
}
