#include "common.h"

#define main library_main
#include "../code_library/bernoulli_numbers.cpp"
#undef main

long long power(long long x, long long n){
    long long res = 1;
    for (x %= MOD; n; n >>= 1, x = x * x % MOD) if (n & 1) res = res * x % MOD;
    return res;
}

/// The classic recurrence sum C(m + 1, j) B_j = 0, independent of the library's Stirling number formula
int main(){
    gen();

    const int n = MAX - 2;
    vector<vector<long long>> binom(n + 2, vector<long long>(n + 2, 0));
    for (int i = 0; i <= n + 1; i++){
        binom[i][0] = 1;
        for (int j = 1; j <= i; j++) binom[i][j] = (binom[i - 1][j - 1] + binom[i - 1][j]) % MOD;
    }

    vector<long long> b(n + 1);
    b[0] = 1;
    for (int m = 1; m <= n; m++){
        long long s = 0;
        for (int j = 0; j < m; j++) s = (s + binom[m + 1][j] * b[j]) % MOD;
        b[m] = (MOD - s) % MOD * power(m + 1, MOD - 2) % MOD;
    }

    for (int i = 0; i <= n; i++) assert(bernoulli[i] == b[i]);

    /// Exact small values: B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, B_12 = -691/2730
    assert(bernoulli[1] == MOD - power(2, MOD - 2));
    assert(bernoulli[2] == power(6, MOD - 2));
    assert(bernoulli[4] == MOD - power(30, MOD - 2));
    assert(bernoulli[12] == (MOD - 691) * power(2730, MOD - 2) % MOD);
    return 0;
}
