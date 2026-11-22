// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_4_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/minimum_weight_cycle.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    MinimumWeightCycle g(n, MinimumWeightCycle::DIRECTED);
    vector<pair<int, int>> edges(m);
    for (auto& [u, v] : edges){
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    WeightCycle c = g.solve();
    if (c.exists){
        int k = c.vertices.size();
        assert(c.weight == k && (int)c.edges.size() == k);
        for (int i = 0; i < k; i++) assert(edges[c.edges[i]] == make_pair(c.vertices[i], c.vertices[(i + 1) % k]));
    }

    printf("%d\n", c.exists ? 1 : 0);
    return 0;
}
