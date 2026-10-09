#include "common.h"

#define main library_main
#include "../code_library/pollard_rho.cpp"
#undef main

/// Independent Miller Rabin with plain __int128 arithmetic, no Montgomery, no sieve
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

long long random_prime(long long lo, long long hi){
    for (;;){
        long long p = stress::rand_int(lo, hi);
        if (reference_prime(p)) return p;
    }
}

/// Factors must multiply back to n, be prime by the reference test and come sorted
void check(long long n, PollardRho& rho, const vector<long long>* expected = nullptr){
    auto f = rho.factorize(n);
    assert(is_sorted(f.begin(), f.end()));
    __int128 prod = 1;
    for (long long p : f){
        assert(reference_prime(p));
        prod *= p;
    }
    assert(prod == n);
    if (expected) assert(f == *expected);
    assert(rho.is_prime(n) == reference_prime(n));
}

int main(){
    PollardRho rho;
    const long long E6 = 1000000, LIM = LLONG_MAX;

    for (long long n = 1; n <= 5000; n++) check(n, rho);
    for (long long n = E6 - 3000; n <= E6 + 3000; n++) check(n, rho);

    for (long long it = 0; it < stress::scaled(5000); it++){
        long long n = it % 2 ? stress::rand_int(1, LIM) : LIM - stress::rand_int(0, 100000);
        check(n, rho);
    }

    /// 299210837 divides the base 1795265022, so that base must be skipped instead of reporting composite
    assert(rho.is_prime(299210837));
    vector<long long> prime_only = {299210837}, with_three = {3, 299210837};
    check(299210837, rho, &prime_only);
    check(3 * 299210837LL, rho, &with_three);

    /// Strong pseudoprimes to small bases and Carmichael numbers
    for (long long n : {561LL, 1105LL, 1729LL, 2047LL, 3215031751LL, 3825123056546413051LL}){
        check(n, rho);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        vector<long long> primes;
        int shape = it % 6;
        if (shape == 0) primes = {random_prime(1000000000, 3037000499), random_prime(1000000000, 3037000499)};
        else if (shape == 1) primes = {random_prime(E6 + 1, 2000000), random_prime(E6 + 1, 2000000), random_prime(E6 + 1, 2000000)};
        else if (shape == 2){
            long long p = random_prime(E6 + 1, 3037000499);
            primes = {p, p};
        }
        else if (shape == 3){
            long long p = random_prime(E6 + 1, 2000000);
            primes = {p, p, random_prime(E6 + 1, 2000000)};
        }
        else if (shape == 4) primes = {random_prime(2, 1000), random_prime(E6 + 1, 1000000000), random_prime(E6 + 1, 1000000000)};
        else primes = {random_prime(2, 50), random_prime(2, 50), random_prime(2, 50), random_prime(1000000000000LL, 1000000000000000LL)};

        __int128 prod = 1;
        for (long long p : primes) prod *= p;
        if (prod > LIM) continue;
        sort(primes.begin(), primes.end());
        check((long long)prod, rho, &primes);
    }

    return 0;
}
