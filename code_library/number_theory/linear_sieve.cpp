/***
 *
 * Linear Sieve
 * Smallest prime factors and primes up to n in O(n), Euler's phi and Mobius mu on request, any multiplicative function
 *
 * Complexity: O(n) time, 4 bytes per entry, 8 with PHI, 5 with MU, 9 with PHI_MU, plus 4 per prime
 *
 * LinearSieve sieve(n): smallest prime factors and primes for 0..n
 * LinearSieve sieve(n, LinearSieve::PHI): also Euler's totient, likewise MU for Mobius mu and PHI_MU for both
 * sieve.factorize(x): prime factors of 1 <= x <= n with multiplicity, ascending, O(log x)
 * sieve.mu(x): Mobius function of 1 <= x <= n, asserts the sieve was built with MU or PHI_MU
 * sieve.multiplicative(f): g[x] for 1 <= x <= n where g is multiplicative and g(p^k) = f(p, k, p^k), g[1] = 1, in O(n)
 *     e.g. divisor count: f = [](long long p, int k, long long pk){ return k + 1; }
 * sieve.phi(x): Euler's totient of 1 <= x <= n, asserts the sieve was built with PHI or PHI_MU
 * sieve.primes: primes <= n in order
 * sieve.spf[x]: smallest prime factor of x >= 2
 *
 * Notes:
 *   - At -O2 n = 1e7 takes about 80 ms, 160 ms with PHI_MU; n = 1e8 needs 0.4 GB, 0.9 GB with PHI_MU
 *   - n = 2^31 - 1 overflows n + 1
 *   - multiplicative() returns a new table of 8 bytes per entry on every call
 *   - For mu alone, mobius_function.cpp needs 1 byte per entry and runs about 1.8x faster
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LinearSieve{
    enum Table{PHI = 1, MU = 2, PHI_MU = 3};

    int n;
    vector<int> spf, primes, phi_table;
    vector<signed char> mu_table;

    LinearSieve(int n, Table tables = Table()) : n(n), spf(n + 1, 0){
        bool want_phi = tables & PHI, want_mu = tables & MU;
        if (want_phi) phi_table.assign(n + 1, 0);
        if (want_mu) mu_table.assign(n + 1, 0);
        if (n >= 1 && want_phi) phi_table[1] = 1;
        if (n >= 1 && want_mu) mu_table[1] = 1;

        for (int i = 2; i <= n; i++){
            if (!spf[i]){
                spf[i] = i, primes.push_back(i);
                if (want_phi) phi_table[i] = i - 1;
                if (want_mu) mu_table[i] = -1;
            }
            for (int p : primes){
                long long x = (long long)i * p;
                if (p > spf[i] || x > n) break;
                spf[x] = p;
                if (want_phi) phi_table[x] = p == spf[i] ? phi_table[i] * p : phi_table[i] * (p - 1);
                if (want_mu) mu_table[x] = p == spf[i] ? 0 : -mu_table[i];
            }
        }
    }

    vector<int> factorize(int x) const{
        assert(1 <= x && x <= n);
        vector<int> res;
        for (; x > 1; x /= spf[x]) res.push_back(spf[x]);
        return res;
    }

    int mu(int x) const{
        assert(1 <= x && x < (int)mu_table.size());  /// mu_table is empty unless the constructor got MU or PHI_MU
        return mu_table[x];
    }

    template <typename F>
    vector<long long> multiplicative(F f) const{
        vector<long long> g(n + 1, 0);

        if (n >= 1) g[1] = 1;
        for (int x = 2; x <= n; x++){
            /// Dividing p out runs 1.2x to 1.7x slower than storing spf exponents, which slowed every sieve by 20%
            int p = spf[x], k = 1, pk = p, rest = x / p;
            while (spf[rest] == p) rest /= p, k++, pk *= p;
            g[x] = g[rest] * f((long long)p, k, (long long)pk);
        }

        return g;
    }

    int phi(int x) const{
        assert(1 <= x && x < (int)phi_table.size());  /// phi_table is empty unless the constructor got PHI or PHI_MU
        return phi_table[x];
    }
};

int main(){
    /// Smallest prime factors, primes and factorizations up to 100
    LinearSieve sieve(100);
    assert(sieve.primes.size() == 25);
    assert(sieve.primes[5] == 13);
    assert(sieve.spf[91] == 7);
    assert(sieve.spf[97] == 97);
    assert((sieve.factorize(60) == vector<int>{2, 2, 3, 5}));
    assert(sieve.factorize(1).empty());

    /// Divisor count and divisor sum from their values on prime powers
    auto divisors = sieve.multiplicative([](long long, int k, long long){ return (long long)k + 1; });
    assert(divisors[12] == 6);
    assert(divisors[64] == 7);

    auto sigma = sieve.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });
    assert(sigma[12] == 28);
    assert(sigma[28] == 56);

    /// Euler's phi and Mobius mu, sieved only when asked for
    LinearSieve both(100, LinearSieve::PHI_MU);
    assert(both.phi(36) == 12);
    assert(both.phi(97) == 96);
    assert(both.mu(30) == -1);
    assert(both.mu(12) == 0);

    return 0;
}
