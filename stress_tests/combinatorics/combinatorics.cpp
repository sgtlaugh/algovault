#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/combinatorics.cpp"
#undef main

bool is_prime(long long p){
    if (p < 2) return false;
    for (long long d = 2; d * d <= p; d++){
        if (p % d == 0) return false;
    }

    return true;
}

long long random_prime(long long lo, long long hi){
    for (;;){
        long long p = stress::rand_int(lo, hi);
        if (is_prime(p)) return p;
    }
}

/// nCr / nPr / inverses against Pascal's triangle and direct products, for primes just above n up to 2^31 - 1
void check_tables(int n, long long mod){
    Combinatorics comb(n, mod);
    vector<long long> row(1, 1);
    for (int i = 0; i <= n; i++){
        for (int k = -1; k <= i + 1; k++){
            long long pascal = (k < 0 || k > i) ? 0 : row[k];
            assert(comb.nCr(i, k) == pascal);

            long long perm = (k < 0 || k > i) ? 0 : 1;
            for (int j = 0; k >= 0 && j < k && k <= i; j++) perm = perm * (i - j) % mod;
            assert(comb.nPr(i, k) == perm);
        }
        if (i >= 1) assert((__int128)comb.inverse(i) * i % mod == 1);
        assert((__int128)comb.factorial(i) * comb.inv_factorial(i) % mod == 1);

        vector<long long> next(i + 2, 1);
        for (int k = 1; k <= i; k++) next[k] = (row[k - 1] + row[k]) % mod;
        row = next;
    }
}

int main(){
    stress::expect_assert("combinatorics/combinatorics.cpp", []{ Combinatorics(10, 1000000007).inverse(11); });
    stress::expect_assert("combinatorics/combinatorics.cpp", []{ Combinatorics(10, 1000000007).inverse(0); });

    const long long INT31 = (1LL << 31) - 1;

    for (int n = 0; n <= 40; n++){
        check_tables(n, random_prime(n + 1, 4 * n + 10));
        check_tables(n, INT31);
    }
    for (long long it = 0; it < stress::scaled(20); it++){
        check_tables(stress::rand_int(50, 400), random_prime(401, INT31));
    }

    return 0;
}
