// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/general_weighted_matching
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/weighted_graph_matching.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    WeightedBlossom wb(n);
    for (int i = 0; i < m; i++){
        int u, v, w;
        if (scanf("%d %d %d", &u, &v, &w) != 3) return 0;
        wb.add_edge(u, v, w);
    }

    long long total = wb.maximum_matching();
    int size = 0;
    for (int v = 0; v < n; v++) size += wb.mate[v] > v;

    printf("%d %lld\n", size, total);
    for (int v = 0; v < n; v++){
        if (wb.mate[v] > v) printf("%d %d\n", v, wb.mate[v]);
    }
    return 0;
}
