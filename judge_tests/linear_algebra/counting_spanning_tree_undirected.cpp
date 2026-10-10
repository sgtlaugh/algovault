// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/counting_spanning_tree_undirected
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/determinant.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    SpanningTrees g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    printf("%lld\n", g.count(998244353));
    return 0;
}
