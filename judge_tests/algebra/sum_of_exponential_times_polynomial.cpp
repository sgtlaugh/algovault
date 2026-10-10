// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sum_of_exponential_times_polynomial
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/lagrange_interpolation.cpp"
#undef main
#define main linear_sieve_main
#include "../../code_library/number_theory/linear_sieve.cpp"
#undef main

const int MOD = 998244353;

long long power(long long b, long long e){
    long long res = 1;
    for (b %= MOD; e; e >>= 1, b = b * b % MOD){
        if (e & 1) res = res * b % MOD;
    }
    return res;
}

int main(){
    long long r, d, n;
    if (scanf("%lld %lld %lld", &r, &d, &n) != 3) return 0;

    if (r == 0){
        puts(n >= 1 && d == 0 ? "1" : "0");
        return 0;
    }

    /// i^d is completely multiplicative, so only primes need a power
    int m = d + 1;
    vector<int> pw(m + 1, 0);
    {
        LinearSieve sieve(m);
        pw[0] = d == 0;
        if (m >= 1) pw[1] = 1;
        for (int x = 2; x <= m; x++){
            int p = sieve.spf[x];
            pw[x] = x == p ? power(x, d) : (long long)pw[x / p] * pw[p] % MOD;
        }
    }

    /// s[k] = sum of r^i i^d over i < k, for k = 0..d+1
    vector<int> s(m + 1, 0);
    for (long long k = 0, rk = 1; k < m; k++, rk = rk * r % MOD) s[k + 1] = (s[k] + rk * pw[k]) % MOD;
    pw.clear(), pw.shrink_to_fit();

    if (r == 1){
        printf("%d\n", Lagrange(s, MOD).interpolate(n));
        return 0;
    }

    /// s[k] = r^k F(k) + c with deg F <= d, the (d+1)-th finite difference of F(k) = (s[k] - c) r^-k vanishes,
    /// so c * (1/r - 1)^(d+1) = sum over k of C(d+1, k) (-1)^(d+1-k) r^-k s[k]
    long long r_inv = power(r, MOD - 2), num = 0, binom = 1, rk = 1;
    vector<int> inv(m + 1, 1);
    for (int k = 2; k <= m; k++) inv[k] = (MOD - (long long)(MOD / k) * inv[MOD % k] % MOD) % MOD;
    for (int k = 0; k <= m; k++){
        long long term = binom * rk % MOD * s[k] % MOD;
        num += (m - k) & 1 ? MOD - term : term;
        if (k < m) binom = binom * (m - k) % MOD * inv[k + 1] % MOD;
        rk = rk * r_inv % MOD;
    }
    long long c = num % MOD * power(power((r_inv + MOD - 1) % MOD, m), MOD - 2) % MOD;

    vector<int> f(m);
    rk = 1;
    for (int k = 0; k < m; k++, rk = rk * r_inv % MOD) f[k] = (s[k] + MOD - c) % MOD * rk % MOD;
    s.clear(), s.shrink_to_fit(), inv.clear(), inv.shrink_to_fit();

    long long rn = power(r, n % (MOD - 1));
    printf("%lld\n", (rn * Lagrange(f, MOD).interpolate(n) + c) % MOD);
    return 0;
}
