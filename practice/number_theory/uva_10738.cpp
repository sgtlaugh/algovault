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
    const int MAXN = 1000000;
    MobiusSieve sieve(MAXN);
    vector<int> mertens(MAXN + 1, 0);
    for (int x = 1; x <= MAXN; x++) mertens[x] = mertens[x - 1] + sieve.mu[x];

    int n;
    while (scanf("%d", &n) == 1 && n) printf("%8d%8d%8d\n", n, sieve.mu[n], mertens[n]);
    return 0;
}
