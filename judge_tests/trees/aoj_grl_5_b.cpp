// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_5_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/rerooting.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    Rerooting<long long> tree(n);
    for (int i = 0; i + 1 < n; i++){
        int s, t, w;
        if (scanf("%d %d %d", &s, &t, &w) != 3) return 0;
        tree.add_edge(s, t, w);
    }

    auto merge = [](long long a, long long b){ return max(a, b); };
    auto apply_edge = [](long long x, int, int, long long w){ return x + w; };
    vector<long long> height = tree.solve(merge, 0LL, apply_edge);
    for (long long h : height) printf("%lld\n", h);
    return 0;
}
