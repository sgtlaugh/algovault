#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/segmented_sieve.cpp"
#undef main

/// Independent deterministic Miller-Rabin with __int128 arithmetic
bool reference_prime(unsigned long long n){
    if (n < 2) return false;
    for (unsigned long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        if (n % p == 0) return n == p;
    }

    auto mulmod = [&](unsigned long long a, unsigned long long b){ return (unsigned long long)((unsigned __int128)a * b % n); };
    auto powmod = [&](unsigned long long a, unsigned long long e){
        unsigned long long r = 1;
        for (; e; e >>= 1, a = mulmod(a, a)){
            if (e & 1) r = mulmod(r, a);
        }
        return r;
    };

    unsigned long long d = n - 1;
    int s = 0;
    for (; !(d & 1); d >>= 1) s++;

    for (unsigned long long a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        unsigned long long x = powmod(a, d);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int i = 1; i < s && composite; i++){
            x = mulmod(x, x);
            if (x == n - 1) composite = false;
        }
        if (composite) return false;
    }

    return true;
}

void check(const SegmentedSieve& sieve, long long L, long long R){
    vector<char> res = sieve.window(L, R);
    assert((long long)res.size() == R - L + 1);
    for (long long x = L; x <= R; x++) assert(res[x - L] == reference_prime(x));
}

/// Windows wider than block_size against an unblocked sieve, too wide for Miller-Rabin on every number
void check_wide(const SegmentedSieve& sieve, const vector<long long>& base, long long L, long long R){
    vector<char> expected(R - L + 1, 1);
    for (long long p : base){
        if (p * p > R) break;
        for (long long j = max(p * p, (L + p - 1) / p * p); j <= R; j += p) expected[j - L] = 0;
    }
    for (long long x = L; x <= min(R, 1LL); x++) expected[x - L] = 0;

    assert(sieve.window(L, R) == expected);
    for (int k = 0; k < 20; k++){
        long long x = stress::rand_int(L, R);
        assert(expected[x - L] == reference_prime(x));
    }
}

int main(){
    const long long MAXR = 100000000000000LL;
    const long long block = SegmentedSieve::block_size;
    SegmentedSieve big(MAXR), near(1003000);

    vector<long long> base;
    vector<char> composite(10000001, 0);
    for (long long p = 2; p <= 10000000; p++){
        if (composite[p]) continue;
        base.push_back(p);
        for (long long j = p * p; j <= 10000000; j += p) composite[j] = 1;
    }

    for (long long L = 0; L <= 60; L++){
        for (long long R = L; R <= 60; R++){
            check(big, L, R);
            check(near, L, R);
        }
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        long long len = stress::rand_int(0, 3000);
        if (it % 3 == 0){
            long long L = stress::rand_int(0, 1000000);
            check(near, L, L + len);
        }
        else{
            long long L = stress::rand_int(0, MAXR - len);
            check(big, L, L + len);
        }
    }

    /// Windows hugging the upper bound and around perfect squares of primes, R = p * p checks the last base prime
    check(big, MAXR - 3000, MAXR);
    for (long long p : {9999991LL, 9999973LL, 1000003LL}){
        check(big, p * p - 1000, p * p + 1000);
        check(SegmentedSieve(p * p), p * p - 1000, p * p);
    }

    /// R - L within 2 of a multiple of block_size catches a last number left to the wrong block
    for (long long k = 1; k <= 3; k++){
        for (long long d = -2; d <= 2; d++){
            long long L = stress::rand_int(0, MAXR - 200000);
            check_wide(big, base, L, L + k * block + d);
            check_wide(near, base, L % 500000, L % 500000 + k * block + d);
        }
    }

    for (long long it = 0; it < stress::scaled(20); it++){
        long long len = stress::rand_int(block, 300000);
        const SegmentedSieve& sieve = it % 4 == 0 ? near : big;
        long long L = it % 3 == 0 ? stress::rand_int(0, 100000) : stress::rand_int(0, sieve.max_r - len);
        check_wide(sieve, base, L, L + len);
    }

    /// The widest windows the header allows, at both ends of the range
    check_wide(big, base, 0, 10000000);
    check_wide(big, base, MAXR - 10000000, MAXR);
    vector<char> first = big.window(0, 10000000);
    assert(count(first.begin(), first.end(), 1) == 664579);

    return 0;
}
