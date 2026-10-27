// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/min_cost_b_flow
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/network_simplex.cpp"
#undef main

void print(__int128 x){
    if (x < 0) putchar('-'), x = -x;
    if (x >= 10) print(x / 10);
    putchar('0' + x % 10);
}

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    NetworkSimplex g(n);
    for (int v = 0; v < n; v++){
        long long b;
        if (scanf("%lld", &b) != 1) return 0;
        g.add_supply(v, b);
    }
    for (int i = 0; i < m; i++){
        int s, t;
        long long low, high, cost;
        if (scanf("%d %d %lld %lld %lld", &s, &t, &low, &high, &cost) != 5) return 0;
        g.add_edge(s, t, low, high, cost);
    }

    auto cost = g.solve();
    if (!cost){
        puts("infeasible");
        return 0;
    }

    print(*cost);
    putchar('\n');
    for (int v = 0; v < n; v++) printf("%lld\n", g.potential[v]);
    for (int i = 0; i < m; i++) printf("%lld\n", g.flow(i));
    return 0;
}
