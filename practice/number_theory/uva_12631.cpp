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

const int MAXP = 100000;
const long long MOD = 1000000007;

long long expo(long long b, long long e){
    long long res = 1;
    for (b %= MOD; e; e >>= 1, b = b * b % MOD){
        if (e & 1) res = res * b % MOD;
    }

    return res;
}

int main(){
    MobiusSieve sieve(MAXP);
    vector<long long> pow2(MAXP + 1, 1), fact(MAXP + 1, 1), inv_fact(MAXP + 1, 1);

    for (int i = 1; i <= MAXP; i++) pow2[i] = pow2[i - 1] * 2 % MOD, fact[i] = fact[i - 1] * i % MOD;
    inv_fact[MAXP] = expo(fact[MAXP], MOD - 2);
    for (int i = MAXP; i > 0; i--) inv_fact[i - 1] = inv_fact[i] * i % MOD;

    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, k;
        long long c, a, b, inc;
        scanf("%d %lld %d %lld %lld %lld", &n, &c, &k, &a, &b, &inc);

        vector<long long> price(n);
        vector<int> freq(MAXP + 1, 0), multiples(MAXP + 1, 0);
        for (int i = 0; i < n; i++){
            price[i] = i == 0 ? a : (i == 1 ? b : 1 + (a * price[i - 2] + b * price[i - 1] + inc) % 100000);
            freq[price[i]]++;
        }
        for (int g = 1; g <= MAXP; g++){
            for (int x = g; x <= MAXP; x += g) multiples[g] += freq[x];
        }

        auto choose = [&](int m){ return m < k ? 0 : fact[m] * inv_fact[k] % MOD * inv_fact[m - k] % MOD; };

        long long any = 0, sized = 0;
        for (int g = 1; g <= MAXP; g++){
            if (!multiples[g]) continue;
            long long all = 0, exact_k = 0;
            for (int j = 1, x = g; x <= MAXP; j++, x += g){
                if (!sieve.mu[j] || !multiples[x]) continue;
                all += sieve.mu[j] * (pow2[multiples[x]] - 1);
                exact_k += sieve.mu[j] * choose(multiples[x]);
            }
            long long ways = (c / g) % MOD;
            any = (any + ways * ((all % MOD + MOD) % MOD)) % MOD;
            sized = (sized + ways * ((exact_k % MOD + MOD) % MOD)) % MOD;
        }
        printf("Case %d: %lld %lld\n", tc, any, sized);
    }

    return 0;
}
