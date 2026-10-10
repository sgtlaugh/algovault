// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_10_D
// competitive-verifier: TLE 2
// competitive-verifier: ERROR 1e-4
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/knuth_optimization.cpp"
#undef main

/// Probabilities have at most 4 decimals, scaling them to exact integers keeps ties exact for the opt monotonicity
long long read_scaled(){
    char buf[32];
    if (scanf("%31s", buf) != 1) return 0;

    long long whole = 0, frac = 0;
    int digits = 0;
    char* p = buf;
    for (; *p && *p != '.'; p++) whole = whole * 10 + (*p - '0');
    if (*p == '.'){
        for (p++; *p && digits < 4; p++, digits++) frac = frac * 10 + (*p - '0');
    }
    for (; digits < 4; digits++) frac *= 10;

    return whole * 10000 + frac;
}

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<long long> p(n + 2, 0), q(n + 1);
    for (int i = 1; i <= n; i++) p[i] = read_scaled();
    for (int i = 0; i <= n; i++) q[i] = read_scaled();

    vector<long long> pre_p(n + 3, 0), pre_q(n + 2, 0);
    for (int i = 0; i < n + 2; i++) pre_p[i + 1] = pre_p[i] + p[i];
    for (int i = 0; i <= n; i++) pre_q[i + 1] = pre_q[i] + q[i];

    /// Dummy key d_i is the leaf [i, i + 1] over points 0..n + 1, splitting at k makes key k the root
    auto dp = knuth_dp<long long>(n + 1, [&](int i, int j){
        return (pre_q[j] - pre_q[i]) + (pre_p[j] - pre_p[i + 1]);
    });

    long long total = dp[0][n + 1] + pre_q[n + 1];
    printf("%.8f\n", total / 10000.0);
    return 0;
}
