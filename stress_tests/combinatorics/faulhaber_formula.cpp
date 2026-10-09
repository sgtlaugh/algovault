#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/faulhaber_formula.cpp"
#undef main

long long power(long long x, long long n){
    long long res = 1;
    for (x %= MOD; n; n >>= 1, x = x * x % MOD) if (n & 1) res = res * x % MOD;
    return res;
}

/// Sum of i^k for i = 1..n by Lagrange interpolation of the degree k + 1 polynomial through x = 0..k + 1
long long lagrange_sum(long long n, int k){
    int m = k + 2;
    vector<long long> ys(m, 0);
    for (int x = 1; x < m; x++) ys[x] = (ys[x - 1] + power(x, k)) % MOD;
    if (n < m) return ys[n];

    long long r = n % MOD;
    vector<long long> pre(m + 1, 1), suf(m + 1, 1), fact(m + 1, 1);
    for (int i = 0; i < m; i++) pre[i + 1] = pre[i] * ((r - i + MOD) % MOD) % MOD;
    for (int i = m - 1; i >= 0; i--) suf[i] = suf[i + 1] * ((r - i + MOD) % MOD) % MOD;
    for (int i = 1; i <= m; i++) fact[i] = fact[i - 1] * i % MOD;

    long long res = 0;
    for (int i = 0; i < m; i++){
        long long term = ys[i] * pre[i] % MOD * suf[i + 1] % MOD * power(fact[i] * fact[m - 1 - i] % MOD, MOD - 2) % MOD;
        res = (res + ((m - 1 - i) % 2 ? MOD - term : term)) % MOD;
    }
    return res;
}

/// Compares against direct summation for small n, Lagrange interpolation for huge n, and the n^k step identity
int main(){
    const int MAXK = 1009;
    Faulhaber f(MAXK);

    for (int max_k = 0; max_k <= 12; max_k++){
        Faulhaber g(max_k);
        for (int k = 0; k <= max_k; k++){
            long long sum = 0;
            for (int n = 0; n <= 60; n++){
                if (n) sum = (sum + power(n, k)) % MOD;
                assert(g.sum(n, k) == sum);
            }
        }
    }

    for (long long it = 0; it < stress::scaled(2000); it++){
        int k = it % 10 ? stress::rand_int(0, 30) : stress::rand_int(0, MAXK);
        int n = stress::rand_int(0, 300);
        long long sum = 0;
        for (int i = 1; i <= n; i++) sum = (sum + power(i, k)) % MOD;
        assert(f.sum(n, k) == sum);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int k = it % 10 ? stress::rand_int(0, 40) : stress::rand_int(0, MAXK);
        if (it == 0) k = MAXK;
        long long n = it % 2 ? stress::rand_int(0, LLONG_MAX) : MOD * stress::rand_int(0, LLONG_MAX / MOD) + stress::rand_int(-3, 3);
        n = max(n, 0LL);
        if (it == 1) n = LLONG_MAX;
        assert(f.sum(n, k) == lagrange_sum(n, k));
    }

    /// Too large to sum: consecutive prefix sums must differ by n^k, including n around multiples of MOD
    for (long long it = 0; it < stress::scaled(20000); it++){
        int k = stress::rand_int(0, it % 10 ? 30 : MAXK);
        long long n = it % 3 ? stress::rand_int(1, LLONG_MAX) : MOD * stress::rand_int(1, LLONG_MAX / MOD - 1) + stress::rand_int(-2, 2);
        assert((f.sum(n, k) - f.sum(n - 1, k) + MOD) % MOD == power(n, k));
    }

    return 0;
}
