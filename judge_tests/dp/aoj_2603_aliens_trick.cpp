// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2603
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/aliens_trick.cpp"
#undef main

int main(){
    int s, n, m;
    if (scanf("%d %d %d", &s, &n, &m) != 3) return 0;

    vector<long long> x(s + 1), a(n);
    for (int i = 1; i <= s; i++){
        if (scanf("%lld", &x[i]) != 1) return 0;
    }
    for (auto& v : a){
        long long t;
        int p;
        if (scanf("%lld %d", &t, &p) != 2) return 0;
        v = t - x[p];
    }

    /// A bus leaving at time D picks up everyone with t - x[p] <= D, so each bus serves a run of the sorted values
    sort(a.begin(), a.end());
    vector<long long> pref(n + 1, 0);
    for (int i = 0; i < n; i++) pref[i + 1] = pref[i] + a[i];

    auto cost = [&](int j, int i){ return a[i - 1] * (i - j) - (pref[i] - pref[j]); };
    auto solver = [&](long long lambda){
        vector<pair<long long, int>> dp(n + 1, {LLONG_MAX, 0});
        dp[0] = {0, 0};
        for (int i = 1; i <= n; i++){
            for (int j = 0; j < i; j++) dp[i] = min(dp[i], {dp[j].first + cost(j, i) + lambda, dp[j].second + 1});
        }

        return dp[n];
    };

    /// The best total is convex and non-increasing in the number of buses, so every slope lies in [0, cost(0, n)]
    printf("%lld\n", aliens_trick(min(n, m), 0LL, cost(0, n), solver));
    return 0;
}
