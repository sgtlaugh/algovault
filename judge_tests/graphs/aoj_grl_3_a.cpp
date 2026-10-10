// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_3_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/articulation_points.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    Graph g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    for (int u : g.get_cuts()) printf("%d\n", u);
    return 0;
}
