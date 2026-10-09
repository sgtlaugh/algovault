#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/segmented_sieve.cpp"
#undef main

/// Independent deterministic Miller Rabin with __int128 arithmetic
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

/// Fresh builds cost O(sqrt(R)) each, so the wrapper is compared on small R and on a few chosen windows
void check(const SegmentedSieve& sieve, long long L, long long R, bool fresh = false){
    vector<char> res = sieve.window(L, R);
    assert((long long)res.size() == R - L + 1);
    long long primes = 0;
    for (long long x = L; x <= R; x++){
        assert(res[x - L] == reference_prime(x));
        primes += res[x - L];
    }
    if (fresh || R <= 10000000000LL) assert(segmented_sieve(L, R) == res);
    if (R <= 10000000000LL) assert(count_primes(L, R) == primes);
}

int main(){
    const long long MAXR = 100000000000000LL;
    SegmentedSieve big(MAXR), near(1003000);
    for (long long L = 0; L <= 60; L++){
        for (long long R = L; R <= 60; R++) check(big, L, R), check(near, L, R);
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
    check(big, MAXR - 3000, MAXR, true);
    for (long long p : {9999991LL, 9999973LL, 1000003LL}){
        check(big, p * p - 1000, p * p + 1000, true);
        check(SegmentedSieve(p * p), p * p - 1000, p * p);
    }
    assert(count_primes(1, 10000000) == 664579);
    return 0;
}
