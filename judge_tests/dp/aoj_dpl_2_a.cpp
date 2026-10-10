// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DPL_2_A
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/hamiltonian_dp.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    HamiltonianGraph g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        long long w;
        if (scanf("%d %d %lld", &u, &v, &w) != 3) return 0;
        g.add_directed_edge(u, v, w);
    }

    long long cost = g.shortest_cycle().first;
    printf("%lld\n", cost == LLONG_MAX ? -1 : cost);
    return 0;
}
