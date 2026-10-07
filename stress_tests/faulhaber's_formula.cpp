#include "common.h"

#define main library_main
#include "../code_library/faulhaber's_formula.c"
#undef main

long long power(long long x, long long n){
    long long res = 1;
    for (x %= MOD; n; n >>= 1, x = x * x % MOD) if (n & 1) res = res * x % MOD;
    return res;
}

int main(){
    generate();

    for (long long it = 0; it < stress::scaled(2000); it++){
        int k = it % 10 ? stress::rand_int(0, 30) : stress::rand_int(0, MAXK - 1);
        int n = stress::rand_int(0, 300);
        long long sum = 0;
        for (int i = 1; i <= n; i++) sum = (sum + power(i, k)) % MOD;
        assert(faulhaber(n, k) == sum);
    }

    /// Too large to sum: consecutive prefix sums must differ by n^k, including n around multiples of MOD
    for (long long it = 0; it < stress::scaled(20000); it++){
        int k = stress::rand_int(0, it % 10 ? 30 : MAXK - 1);
        long long n = it % 3 ? stress::rand_int(1, 1000000000000000000LL) : MOD * stress::rand_int(1, 1000000) + stress::rand_int(-2, 2);
        assert((faulhaber(n, k) - faulhaber(n - 1, k) + MOD) % MOD == power(n, k));
    }
    return 0;
}
