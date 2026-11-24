// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2603
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/convex_hull_trick.cpp"
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

    /// dp[i] = min over j of (prev[j] + pref[j]) - j * a[i - 1], plus i * a[i - 1] - pref[i]: line j, slope -j
    const long long INF = LLONG_MAX;
    vector<long long> prev(n + 1, INF), dp(n + 1);
    prev[0] = 0;
    for (int k = 0; k < min(n, m); k++){
        ConvexHullTrick<long long> hull;
        bool any = false;
        dp[0] = INF;
        for (int i = 1; i <= n; i++){
            if (prev[i - 1] != INF) hull.add_line(-(i - 1), prev[i - 1] + pref[i - 1]), any = true;
            dp[i] = any ? hull.query_increasing(a[i - 1]) + i * a[i - 1] - pref[i] : INF;
        }
        swap(prev, dp);
    }

    printf("%lld\n", prev[n]);
    return 0;
}
