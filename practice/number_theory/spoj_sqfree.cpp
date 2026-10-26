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
 *   - For spf, phi or any multiplicative function use linear_sieve.cpp, at about 8x the memory and 3x the time
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

#define MAX 10000001

int mertens[MAX];

int main(){
    int t, i;
    long long n, d, e, q, res;

    MobiusSieve sieve(MAX);
    for (i = 1; i < MAX; i++) mertens[i] = mertens[i - 1] + sieve.mu[i];

    scanf("%d", &t);
    while (t--){
        scanf("%lld", &n);
        res = 0;
        for (d = 1; (d * d) <= n; d = e + 1){
            q = n / (d * d);
            e = sqrt((double)(n / q));
            res += q * (mertens[e] - mertens[d - 1]);
        }
        printf("%lld\n", res);
    }
    return 0;
}
