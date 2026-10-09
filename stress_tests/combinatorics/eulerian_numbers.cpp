#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/eulerian_numbers.c"
#undef main

long long power(long long x, long long n){
    long long res = 1;
    for (x %= MOD; n; n >>= 1, x = x * x % MOD) if (n & 1) res = res * x % MOD;
    return res;
}

/// A(n, k) = sum (-1)^j C(n + 1, j) (k + 1 - j)^n over j <= k, checked for the whole row n
void check_row(int n){
    vector<long long> binom(n + 2), pw(n + 2);
    binom[0] = 1;
    for (int j = 1; j <= n + 1; j++) binom[j] = binom[j - 1] * (n + 2 - j) % MOD * power(j, MOD - 2) % MOD;
    for (int x = 0; x <= n + 1; x++) pw[x] = power(x, n);

    for (int k = 0; k < n; k++){
        long long a = 0;
        for (int j = 0; j <= k; j++) a = (a + (j & 1 ? MOD - 1 : 1) * binom[j] % MOD * pw[k + 1 - j]) % MOD;
        assert(dp[n][k] == a);
    }
    assert(dp[n][n] == (n == 0));
}

int main(){
    generate();

    /// Ascents counted over every permutation
    for (int n = 1; n <= 8; n++){
        vector<int> p(n), count(n + 1, 0);
        iota(p.begin(), p.end(), 1);
        do {
            int ascents = 0;
            for (int i = 1; i < n; i++) ascents += p[i] > p[i - 1];
            count[ascents]++;
        } while (next_permutation(p.begin(), p.end()));
        for (int k = 0; k <= n; k++) assert(dp[n][k] == count[k]);
    }

    for (int n = 0; n <= 60; n++) check_row(n);
    for (long long it = 0; it < stress::scaled(6); it++) check_row(stress::rand_int(it % 2 ? 61 : MAX - 100, MAX - 1));

    return 0;
}
