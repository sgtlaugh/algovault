// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_3_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/bridge.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    static Graph g(n);  /// a temporary Graph is several MB, on the stack it leaves too little room for the recursive DFS
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    vector<pair<int, int>> bridges;
    for (auto& b : g.get_bridges()) bridges.push_back({min(b.u, b.v), max(b.u, b.v)});
    sort(bridges.begin(), bridges.end());

    for (auto [u, v] : bridges) printf("%d %d\n", u, v);
    return 0;
}
