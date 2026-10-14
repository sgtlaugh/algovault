/***
 *
 * Linear Sieve
 * Smallest prime factors, primes, Euler's phi, Mobius mu and any multiplicative function for 1..n in one pass
 *
 * Complexity: O(n) time, O(n) memory, about 20 bytes per entry (200 MB at n = 1e7)
 * For mu alone, mobius_function.cpp needs 1 byte per entry and runs about 5x faster
 *
 * LinearSieve sieve(n): tables for 0..n
 * sieve.spf[x]: smallest prime factor of x >= 2, sieve.primes: primes <= n in order
 * sieve.phi[x], sieve.mu[x]: Euler's totient and Mobius function for x >= 1
 * sieve.factorize(x): prime factors of 1 <= x <= n with multiplicity, ascending, O(log x)
 * sieve.multiplicative(f): g[x] for 1 <= x <= n where g is multiplicative and g(p^k) = f(p, k, p^k), g[1] = 1
 *     e.g. divisor count: f = [](long long p, int k, long long pk){ return k + 1; }
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LinearSieve{
    int n;
    vector<int> spf, primes, phi, mu, power, exponent;  /// power[x] = exact power of spf[x] dividing x, exponent[x] = its exponent

    LinearSieve(int n) : n(n), spf(n + 1, 0), phi(n + 1, 0), mu(n + 1, 0), power(n + 1, 0), exponent(n + 1, 0){
        if (n >= 1) phi[1] = mu[1] = 1;

        for (int i = 2; i <= n; i++){
            if (!spf[i]){
                spf[i] = i, primes.push_back(i);
                phi[i] = i - 1, mu[i] = -1, power[i] = i, exponent[i] = 1;
            }
            for (int p : primes){
                long long x = (long long)i * p;
                if (p > spf[i] || x > n) break;
                spf[x] = p;
                if (p == spf[i]){
                    phi[x] = phi[i] * p, mu[x] = 0;
                    power[x] = power[i] * p, exponent[x] = exponent[i] + 1;
                }
                else{
                    phi[x] = phi[i] * (p - 1), mu[x] = -mu[i];
                    power[x] = p, exponent[x] = 1;
                }
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
            int rest = x / power[x];
            g[x] = g[rest] * f((long long)spf[x], exponent[x], (long long)power[x]);
        }
        return g;
    }
};

const int MAXD = 200000;

int main(){
    LinearSieve sieve(MAXD);
    /// count[d] = fractions printed for denominators up to d, 0/1 included
    vector<long long> count(MAXD + 1, 1);
    for (int d = 1; d <= MAXD; d++) count[d] = count[d - 1] + sieve.phi[d];

    long long k;
    while (scanf("%lld", &k) == 1 && k){
        if (k == 1){
            puts("0/1");
            continue;
        }

        int d = lower_bound(count.begin(), count.end(), k) - count.begin();
        long long rank = k - count[d - 1];

        /// Numerators up to x coprime to d: sum of mu(q) * (x / q) over the squarefree divisors q of d
        vector<int> divisors = {1}, factors = sieve.factorize(d);
        factors.erase(unique(factors.begin(), factors.end()), factors.end());
        for (int p : factors){
            for (int i = (int)divisors.size() - 1; i >= 0; i--) divisors.push_back(divisors[i] * p);
        }
        auto coprime_upto = [&](int x){
            int res = 0;
            for (int q : divisors) res += sieve.mu[q] * (x / q);
            return res;
        };

        int lo = 1, hi = d;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (coprime_upto(mid) >= rank) hi = mid;
            else lo = mid + 1;
        }
        printf("%d/%d\n", lo, d);
    }
    return 0;
}
