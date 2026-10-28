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

bool suspect(unsigned long long b, unsigned long long n){
    unsigned long long u = n - 1, x = 1;
    int t = 0;

    while (!(u & 1)) u >>= 1, t++;
    for (b %= n; u; u >>= 1, b = b * b % n){
        if (u & 1) x = x * b % n;
    }
    while (t--){
        unsigned long long y = x * x % n;
        if (y == 1 && x != 1 && x != n - 1) return false;
        x = y;
    }

    return x == 1;
}

int main(){
    SegmentedSieve sieve(4294967295LL);
    long long base, lo, hi;

    for (int tc = 0; scanf("%lld %lld %lld", &base, &lo, &hi) == 3 && base; tc++){
        vector<char> prime = sieve.window(lo, hi);
        vector<long long> fails;
        int composites = 0;
        for (long long n = lo | 1; n <= hi; n += 2){
            if (prime[n - lo]) continue;
            composites++;
            if (suspect(base, n)) fails.push_back(n);
        }

        if (tc) puts("");
        printf("There are %d odd non-prime numbers between %lld and %lld.\n", composites, lo, hi);
        if (fails.empty()) printf("There are no failures in base %lld.\n", base);
        else{
            printf("%d suspects fail in base %lld:\n", (int)fails.size(), base);
            for (long long n : fails) printf("%lld\n", n);
        }
    }

    return 0;
}
