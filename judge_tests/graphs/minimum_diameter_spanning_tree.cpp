// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/minimum_diameter_spanning_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/absolute_center.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    AbsoluteCenter g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        long long w;
        if (scanf("%d %d %lld", &u, &v, &w) != 3) return 0;
        g.add_edge(u, v, w);
    }

    AbsoluteCenterResult r = g.solve();
    printf("%lld\n", r.diameter);
    for (int id : r.tree) printf("%d ", id);
    printf("\n");
    return 0;
}
