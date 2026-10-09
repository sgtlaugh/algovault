/***
 *
 * Segmented Sieve
 * Primality of every number in a window [L, R] far from zero
 *
 * Complexity: O(sqrt(R) log log R) once for the base primes, then O((R - L + 1) log log R + pi(sqrt(R))) per window
 * O(sqrt(R) + R - L) memory
 *
 * SegmentedSieve sieve(max_r): base primes up to sqrt(max_r), for max_r <= 1e14
 *   sieve.window(L, R): res[i] = 1 if L + i is prime, for 0 <= L <= R <= max_r and R - L <= 1e7
 *   Build it once for many windows, at R around 1e14 the base primes dominate a single window
 * segmented_sieve(L, R): one window, builds a fresh SegmentedSieve(R)
 * count_primes(L, R), primes_in_range(L, R): convenience wrappers over segmented_sieve
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SegmentedSieve{
    long long max_r;
    vector<int> primes;

    /// sqrtl is exact on integers up to 1e14, so the floor needs no correction
    SegmentedSieve(long long max_r) : max_r(max_r){
        assert(0 <= max_r && max_r <= 100000000000000LL);
        int limit = sqrtl((long double)max_r);
        vector<char> composite(limit + 1, 0);
        for (int p = 2; p <= limit; p++){
            if (composite[p]) continue;
            primes.push_back(p);
            for (long long j = (long long)p * p; j <= limit; j += p) composite[j] = 1;
        }
    }

    vector<char> window(long long L, long long R) const{
        assert(0 <= L && L <= R && R <= max_r && R - L <= 10000000);
        vector<char> res(R - L + 1, 1);
        for (long long p : primes){
            if (p * p > R) break;
            for (long long j = max(p * p, (L + p - 1) / p * p); j <= R; j += p) res[j - L] = 0;
        }
        for (long long x = L; x <= min(R, 1LL); x++) res[x - L] = 0;
        return res;
    }
};

vector<char> segmented_sieve(long long L, long long R){
    return SegmentedSieve(R).window(L, R);
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

    SegmentedSieve sieve(1000000);
    vector<char> window = sieve.window(999980, 1000000);
    assert(window[999983 - 999980] && count(window.begin(), window.end(), 1) == 1);
    window = sieve.window(0, 30);
    assert(count(window.begin(), window.end(), 1) == 10);
    return 0;
}
