#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/modular_sqrt.cpp"
#undef main

bool brute_is_prime(long long n){
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; d++){
        if (n % d == 0) return false;
    }

    return true;
}

/// Deterministic Miller-Rabin with the first 12 prime bases (exact below 3.3e24) on plain __int128 arithmetic
bool is_prime(long long n){
    auto power = [](__int128 x, long long e, long long m){
        __int128 res = 1;
        for (x %= m; e; e >>= 1, x = x * x % m) if (e & 1) res = res * x % m;
        return (long long)res;
    };

    if (n < 2) return false;
    for (long long q : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        if (n % q == 0) return n == q;
    }

    long long d = n - 1;
    int s = 0;
    for (; d % 2 == 0; d /= 2) s++;
    for (long long base : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
        __int128 x = power(base, d, n);
        for (int i = 1; i < s && x != 1 && x != n - 1; i++) x = x * x % n;
        if (x != 1 && x != n - 1) return false;
    }

    return true;
}

long long random_prime(long long lo, long long hi){
    while (1){
        long long n = stress::rand_int(lo, hi);
        if (is_prime(n)) return n;
    }
}

/// Square roots against a scan over every residue for small primes, then defining properties on random 62-bit primes
int main(){
    for (long long p = 2; p < 1000; p++){
        if (!brute_is_prime(p)) continue;
        vector<long long> smallest(p, -1);
        for (long long x = p - 1; x >= 0; x--) smallest[x * x % p] = x;
        for (long long a = -p; a < 2 * p; a++) assert(sqrt_mod(a, p) == smallest[(a % p + p) % p]);
    }

    /// Random 62-bit primes, the largest prime below 2^62, and two with a deep 2-adic p - 1 for Tonelli-Shanks
    const vector<long long> fixed_primes = {998244353, 29 * (1LL << 57) + 1, (1LL << 62) - 57};
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long p = it % 3 ? random_prime(1LL << 61, (1LL << 62) - 1) : fixed_primes[it / 3 % 3];
        long long x = stress::rand_int(1, p - 1), a = (unsigned __int128)x * x % p;
        assert(sqrt_mod(a, p) == min(x, p - x));
        assert(sqrt_mod(a - p, p) == min(x, p - x));

        long long n = stress::rand_int(1, p - 1);
        long long r = sqrt_mod(n, p);
        if (r != -1) assert((unsigned __int128)r * r % p == (unsigned long long)n && 2 * r < p);
        else assert(sqrt_mod((unsigned __int128)n * a % p, p) == -1);
    }

    /// a = x^(2^40) mod 29 * 2^57 + 1 makes Tonelli-Shanks shift by up to 2^39, random squares take about 2^32 tries
    const long long Q = 29 * (1LL << 57) + 1;
    for (long long it = 0; it < stress::scaled(300); it++){
        long long x = stress::rand_int(1, Q - 1), root = x;
        for (int i = 0; i < 39; i++) root = (unsigned __int128)root * root % Q;
        long long a = (unsigned __int128)root * root % Q;
        assert(sqrt_mod(a, Q) == min(root, Q - root));
    }

    return 0;
}
