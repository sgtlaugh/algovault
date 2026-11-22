// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2603
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/divide_conquer_dp.cpp"
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
    printf("%lld\n", divide_conquer_dp<long long>(n, min(n, m), cost)[n]);
    return 0;
}
