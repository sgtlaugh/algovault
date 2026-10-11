// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_6_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/mcmf_spfa.cpp"
#undef main

int main(){
    int n, m, f;
    if (scanf("%d %d %d", &n, &m, &f) != 3) return 0;

    MCMFSpfa g(n);
    for (int i = 0; i < m; i++){
        int u, v, c, d;
        if (scanf("%d %d %d %d", &u, &v, &c, &d) != 4) return 0;
        g.add_edge(u, v, c, d);
    }

    auto [flow, cost] = g.solve(0, n - 1, f);
    printf("%lld\n", flow == f ? cost : -1LL);
    return 0;
}
