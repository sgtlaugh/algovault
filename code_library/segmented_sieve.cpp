/***
 *
 * Segmented Sieve
 * Primality of every number in a window [L, R] far from zero
 *
 * Complexity: O(sqrt(R) log log R + (R - L + 1) log log R), O(sqrt(R) + R - L) memory
 *
 * segmented_sieve(L, R): res[i] = 1 if L + i is prime, for 0 <= L <= R <= 1e14 and R - L <= 1e7
 * count_primes(L, R), primes_in_range(L, R): convenience wrappers
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

vector<char> segmented_sieve(long long L, long long R){
    assert(0 <= L && L <= R && R <= 100000000000000LL && R - L <= 10000000);

    long long limit = sqrtl((long double)R);
    while (limit * limit > R) limit--;
    while ((limit + 1) * (limit + 1) <= R) limit++;

    vector<char> small(limit + 1, 1), res(R - L + 1, 1);
    for (long long p = 2; p <= limit; p++){
        if (!small[p]) continue;
        for (long long j = p * p; j <= limit; j += p) small[j] = 0;
        for (long long j = max(p * p, (L + p - 1) / p * p); j <= R; j += p) res[j - L] = 0;
    }
    for (long long x = L; x <= min(R, 1LL); x++) res[x - L] = 0;
    return res;
}

long long count_primes(long long L, long long R){
    vector<char> sieve = segmented_sieve(L, R);
    return count(sieve.begin(), sieve.end(), 1);
}

vector<long long> primes_in_range(long long L, long long R){
    vector<char> sieve = segmented_sieve(L, R);
    vector<long long> res;
    for (long long i = 0; i <= R - L; i++){
        if (sieve[i]) res.push_back(L + i);
    }
    return res;
}

int main(){
    assert((primes_in_range(0, 30) == vector<long long>{2, 3, 5, 7, 11, 13, 17, 19, 23, 29}));
    assert((primes_in_range(1, 2) == vector<long long>{2}));
    assert(primes_in_range(0, 1).empty());
    assert(primes_in_range(24, 28).empty());
    assert((primes_in_range(1000000000, 1000000021) == vector<long long>{1000000007, 1000000009, 1000000021}));
    assert(count_primes(1, 100) == 25);
    assert(count_primes(1, 1000000) == 78498);
    assert(count_primes(99999999999974LL, 100000000000000LL) == 0);
    assert(primes_in_range(99999999999970LL, 100000000000000LL).back() == 99999999999973LL);
    return 0;
}
