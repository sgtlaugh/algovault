/***
 *
 * Mobius Function
 * Mobius mu for 1..n, mu[1] = 1, mu[x] = 0 if a squared prime divides x, else (-1)^(number of prime factors)
 *
 * Complexity: O(n) time, n + 1 bytes for mu plus 4 bytes per prime
 *
 * MobiusSieve sieve(n): sieve.mu[x] for 1 <= x <= n, sieve.primes: primes <= n in order
 * mobius_by_divisors(n): the same mu table in O(n log n) from sum over d | x of mu[d] = [x == 1], less code
 *
 * linear_sieve.cpp also gives mu, along with spf, phi and any multiplicative function, at about 20 bytes per entry
 * When only mu is needed this uses about 12x less memory and runs about 5x faster, n = 1e8 fits in about 135 MB
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MobiusSieve{
    vector<signed char> mu;
    vector<int> primes;

    MobiusSieve(int n) : mu(n + 1, 2){  /// 2 marks numbers no smaller prime has reached, which are exactly the primes
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

vector<signed char> mobius_by_divisors(int n){
    vector<signed char> mu(n + 1, 0);
    if (n >= 1) mu[1] = 1;

    for (int i = 1; i <= n; i++){
        for (int j = 2 * i; j <= n; j += i) mu[j] -= mu[i];
    }
    return mu;
}

int main(){
    MobiusSieve sieve(10);
    assert(vector<int>(sieve.mu.begin(), sieve.mu.end()) == vector<int>({0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1}));
    assert(sieve.primes == vector<int>({2, 3, 5, 7}));
    assert(mobius_by_divisors(10) == sieve.mu);

    /// Mertens function M(10^6) = 212
    MobiusSieve big(1000000);
    assert(accumulate(big.mu.begin(), big.mu.end(), 0) == 212);
    assert(big.primes.size() == 78498);

    MobiusSieve empty(0);
    assert(empty.mu.size() == 1 && empty.primes.empty());
    return 0;
}
