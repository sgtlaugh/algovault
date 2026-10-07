#include "common.h"

#define main library_main
#include "../code_library/miller_rabin.cpp"
#undef main

/// Independent reference: Miller-Rabin with the first 12 prime bases is exact below 3.3 * 10^24
bool reference_is_prime(unsigned long long n){
    if (n < 2) return false;
    const int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int p : bases) if (n % p == 0) return n == (unsigned long long)p;

    unsigned long long d = n - 1;
    int s = 0;
    while (!(d & 1)) d >>= 1, s++;
    auto mul = [&](unsigned long long a, unsigned long long b){ return (unsigned long long)((unsigned __int128)a * b % n); };

    for (int a : bases){
        unsigned long long x = 1, b = a;
        for (unsigned long long e = d; e; e >>= 1, b = mul(b, b)) if (e & 1) x = mul(x, b);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s && composite; r++) composite = (x = mul(x, x)) != n - 1;
        if (composite) return false;
    }
    return true;
}

int main(){
    const int N = 200000;
    vector<bool> sieve(N, true);
    sieve[0] = sieve[1] = false;
    for (int i = 2; i * i < N; i++) if (sieve[i]) for (int j = i * i; j < N; j += i) sieve[j] = false;
    for (int n = -5; n < N; n++) assert(prm::is_prime(n) == (n >= 0 && sieve[n]));

    const long long special[] = {2047, 1373653, 25326001, 3215031751LL, 2152302898747LL, 3474749660383LL, 341550071728321LL,
                                 3825123056546413051LL, 561, 1105, 1729, 299210837, 407521, 2147483647, 2147483649LL,
                                 4294967291LL, 4294967295LL, 4294967311LL, 9223372036854775783LL, LLONG_MAX};
    for (long long n : special) assert(prm::is_prime(n) == reference_is_prime(n));

    for (long long it = 0; it < stress::scaled(1000000); it++){
        long long n = (long long)(stress::rng()() >> stress::rand_int(1, 62));
        assert(prm::is_prime(n) == reference_is_prime(n));
    }

    /// Semiprimes of two large factors are the hardest composites for Miller-Rabin
    for (long long it = 0; it < stress::scaled(100000); it++){
        long long p = stress::rand_int(2, 3000000000LL), q = stress::rand_int(2, 3000000000LL);
        if ((unsigned __int128)p * q > LLONG_MAX) continue;
        assert(!prm::is_prime(p * q));
    }
    return 0;
}
