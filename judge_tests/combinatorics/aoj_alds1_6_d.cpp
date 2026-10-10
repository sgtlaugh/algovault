// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_6_D
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/permutation_cycles.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<int> w(n);
    for (auto& x: w){
        if (scanf("%d", &x) != 1) return 0;
    }

    vector<int> sorted_w = w, p(n);
    sort(sorted_w.begin(), sorted_w.end());
    for (int i = 0; i < n; i++) p[i] = lower_bound(sorted_w.begin(), sorted_w.end(), w[i]) - sorted_w.begin();

    /// fix a cycle with its own minimum, or swap the global minimum in and back out
    long long res = 0, global_min = sorted_w[0];
    for (const auto& cycle: permutation_cycles(p)){
        long long sum = 0, cycle_min = LLONG_MAX, len = cycle.size();
        for (int i: cycle) sum += w[i], cycle_min = min(cycle_min, (long long)w[i]);
        res += min(sum + (len - 2) * cycle_min, sum + cycle_min + (len + 1) * global_min);
    }

    printf("%lld\n", res);
    return 0;
}
