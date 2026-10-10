// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_6_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/mcmf_spfa.cpp"
#undef main

int main(){
    int n, m, f;
    if (scanf("%d %d %d", &n, &m, &f) != 3) return 0;

    /// FlowGraph(n) covers nodes 0..n, so node n is free for a super source that limits the flow to f
    static FlowGraph<long long> g(n);
    g.add_edge(n, 0, f, 0);
    for (int i = 0; i < m; i++){
        int u, v, c, d;
        if (scanf("%d %d %d %d", &u, &v, &c, &d) != 4) return 0;
        g.add_edge(u, v, c, d);
    }

    auto [cost, flow] = g.mincost_maxflow(n, n - 1);
    printf("%lld\n", flow == f ? cost : -1LL);
    return 0;
}
