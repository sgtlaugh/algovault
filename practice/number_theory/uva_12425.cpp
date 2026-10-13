/***
 *
 * Divisors of a single number
 *
 * Complexity: O(sqrt(n)) to factorize by trial division, O(d(n) log d(n)) to list and sort the divisors
 *
 * divisors(n): every divisor of 1 <= n <= 1e14 in ascending order
 * divisors_from_factors(f): the same from prime factors with multiplicity, e.g. pollard_rho.cpp's factorize(n)
 *     use this for n up to 2^63, products stay below n so nothing overflows
 *
 * Tables of the divisor count or divisor sum for every n up to a limit: linear_sieve.cpp's multiplicative()
 * Every divisor of every n up to a limit: all_divisors.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<long long> divisors_from_factors(vector<long long> factors){
    sort(factors.begin(), factors.end());
    vector<long long> res = {1};

    for (size_t i = 0; i < factors.size(); ){
        size_t j = i;
        while (j < factors.size() && factors[j] == factors[i]) j++;
        size_t old = res.size();
        long long power = 1;
        for (size_t k = i; k < j; k++){
            power *= factors[i];
            for (size_t t = 0; t < old; t++) res.push_back(res[t] * power);
        }
        i = j;
    }

    sort(res.begin(), res.end());
    return res;
}

vector<long long> divisors(long long n){
    assert(1 <= n && n <= 100000000000000LL);
    vector<long long> factors;

    for (long long p = 2; p * p <= n; p++){
        for (; n % p == 0; n /= p) factors.push_back(p);
    }
    if (n > 1) factors.push_back(n);

    return divisors_from_factors(factors);
}

int main(){
    int t, q;
    long long n, x;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        scanf("%lld %d", &n, &q);
        vector<long long> div = divisors(n), primes;
        for (long long d : div){
            if (d > 1 && all_of(primes.begin(), primes.end(), [&](long long p){ return d % p; })) primes.push_back(d);
        }

        /// Exactly phi(n / d) numbers I <= n have gcd(I, n) = d, so prefix[i] counts the gcds up to div[i]
        int k = div.size();
        vector<long long> prefix(k + 1, 0);
        for (int i = 0; i < k; i++){
            long long m = n / div[i], phi = m;
            for (long long p : primes){
                if (m % p == 0) phi = phi / p * (p - 1);
            }
            prefix[i + 1] = prefix[i] + phi;
        }

        printf("Case %d\n", tc);
        while (q--){
            scanf("%lld", &x);
            printf("%lld\n", prefix[upper_bound(div.begin(), div.end(), x) - div.begin()]);
        }
    }
    return 0;
}
