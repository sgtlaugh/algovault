// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2603
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/concave_1d1d_dp.cpp"
#undef main

#define main aliens_main
#include "../../code_library/dp/aliens_trick.cpp"
#undef main

int main(){
    int s, n, m;
    if (scanf("%d %d %d", &s, &n, &m) != 3) return 0;

    vector<long long> x(s + 1), a(n);
    for (int i = 1; i <= s; i++) scanf("%lld", &x[i]);
    for (auto& v : a){
        long long t;
        int p;
        scanf("%lld %d", &t, &p);
        v = t - x[p];
    }

    /// A bus leaving at time D picks up everyone with t - x[p] <= D, so each bus serves a run of the sorted values
    sort(a.begin(), a.end());
    vector<long long> pref(n + 1, 0);
    for (int i = 0; i < n; i++) pref[i + 1] = pref[i] + a[i];

    auto cost = [&](int j, int i){ return a[i - 1] * (i - j) - (pref[i] - pref[j]); };
    auto solver = [&](long long lambda){
        Concave1D1D<long long> res = concave_1d1d_dp_opt(n, 0LL, [&](int i, int x){ return cost(i, x) + lambda; });
        int buses = 0;
        for (int v = n; v > 0; v = res.opt[v]) buses++;
        return make_pair(res.dp[n], buses);
    };

    /// The best total is convex and non-increasing in the number of buses, so every slope lies in [0, cost(0, n)]
    printf("%lld\n", aliens_trick(min(n, m), 0LL, cost(0, n), solver));
    return 0;
}
