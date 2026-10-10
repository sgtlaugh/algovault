// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_6_A
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/maxflow.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    FlowGraph sparse(n, 0, n - 1);
    DenseFlowGraph dense(n, 0, n - 1);
    for (int i = 0; i < m; i++){
        int u, v;
        long long c;
        if (scanf("%d %d %lld", &u, &v, &c) != 3) return 0;
        sparse.add_directed_edge(u, v, c);
        dense.add_directed_edge(u, v, c);
    }

    long long flow = sparse.maxflow();
    assert(dense.maxflow() == flow);
    printf("%lld\n", flow);
    return 0;
}
