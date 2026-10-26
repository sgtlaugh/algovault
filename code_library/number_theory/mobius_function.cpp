/***
 *
 * Mobius Function
 * Mobius mu for 1..n, mu[1] = 1, mu[x] = 0 if a squared prime divides x, else (-1)^(number of prime factors)
 *
 * Complexity: O(n) time, n + 1 bytes for mu plus 4 bytes per prime
 *
 * MobiusSieve sieve(n): sieves 0..n for n >= 0
 * sieve.mu[x]: Mobius function of 1 <= x <= n
 * sieve.primes: primes <= n in order
 *
 * Notes:
 *   - mu holds signed char, so cast it to int before printing with cout
 *   - For sum mu(d) g(n / d), take d <= s = sqrt(n) singly, then d in (n / (q + 1), n / q] for each q <= n / (s + 1)
 *   - That loop divides once per block, about 4x faster at n = 1e7 than the usual blocks ending at n / (n / l)
 *   - Inverting g(x) = sum of f(d) over d | x needs no mu: divisor_mobius in gcd_lcm_convolution.cpp, O(n log log n)
 *   - For spf, phi or any multiplicative function use linear_sieve.cpp: 4 to 9 bytes per entry, 1.3x to 3x the time
 *   - n = 1e8 fits in about 135 MB
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MobiusSieve{
    vector<signed char> mu;
    vector<int> primes;

    MobiusSieve(int n) : mu(n + 1, 2){  /// 2 marks numbers no smaller prime has reached, which are exactly the primes
        assert(n >= 0);
        mu[0] = 0;
        if (n >= 1) mu[1] = 1;

        for (int i = 2; i <= n; i++){
            if (mu[i] == 2) mu[i] = -1, primes.push_back(i);
            for (int p : primes){
                long long x = (long long)i * p;
                if (x > n) break;
                if (i % p == 0){
                    mu[x] = 0;
                    break;
                }
                mu[x] = -mu[i];
            }
        }
    }
};

int main(){
    /// mu and the primes up to 10
    MobiusSieve sieve(10);
    assert(vector<int>(sieve.mu.begin(), sieve.mu.end()) == vector<int>({0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1}));
    assert(sieve.primes == vector<int>({2, 3, 5, 7}));

    /// Mertens function M(10^6) = 212
    MobiusSieve big(1000000);
    assert(accumulate(big.mu.begin(), big.mu.end(), 0) == 212);
    assert(big.primes.size() == 78498);

    return 0;
}
