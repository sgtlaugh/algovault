/***
 *
 * Linear Sieve
 * Smallest prime factors, primes, Euler's phi, Mobius mu and any multiplicative function for 1..n in one pass
 *
 * Complexity: O(n) time, O(n) memory, about 20 bytes per entry (200 MB at n = 1e7)
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

int main(){
    LinearSieve sieve(100);
    assert((vector<int>(sieve.primes.begin(), sieve.primes.begin() + 6) == vector<int>{2, 3, 5, 7, 11, 13}));
    assert(sieve.primes.size() == 25);
    assert(sieve.spf[91] == 7 && sieve.spf[97] == 97 && sieve.spf[64] == 2);

    assert(sieve.phi[1] == 1 && sieve.phi[8] == 4 && sieve.phi[9] == 6 && sieve.phi[36] == 12 && sieve.phi[97] == 96);
    assert(sieve.mu[1] == 1 && sieve.mu[6] == 1 && sieve.mu[30] == -1 && sieve.mu[12] == 0 && sieve.mu[97] == -1);
    assert((sieve.factorize(60) == vector<int>{2, 2, 3, 5}) && sieve.factorize(1).empty());

    auto divisors = sieve.multiplicative([](long long, int k, long long){ return (long long)k + 1; });
    assert(divisors[1] == 1 && divisors[12] == 6 && divisors[36] == 9 && divisors[97] == 2 && divisors[64] == 7);

    auto sigma = sieve.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });
    assert(sigma[6] == 12 && sigma[12] == 28 && sigma[28] == 56);

    LinearSieve tiny(1);
    assert(tiny.primes.empty() && tiny.phi[1] == 1);

    return 0;
}
