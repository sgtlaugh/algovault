// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/bipartite_edge_coloring
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/edge_coloring.cpp"
#undef main

int main(){
    int l, r, m;
    if (scanf("%d %d %d", &l, &r, &m) != 3) return 0;

    BipartiteEdgeColoring g(l, r);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    auto color = g.solve();
    int k = *max_element(color.begin(), color.end()) + 1;

    string out = to_string(k) + "\n";
    for (int c : color) out += to_string(c) + "\n";
    fputs(out.c_str(), stdout);
    return 0;
}
