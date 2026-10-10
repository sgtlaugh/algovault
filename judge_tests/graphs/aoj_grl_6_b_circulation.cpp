// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_6_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/min_cost_circulation.cpp"
#undef main

int main(){
    int n, m;
    long long f;
    if (scanf("%d %d %lld", &n, &m, &f) != 3) return 0;
    MinCostCirculation g(n);
    long long big = 1;
    for (int i = 0; i < m; i++){
        int u, v;
        long long c, d;
        if (scanf("%d %d %lld %lld", &u, &v, &c, &d) != 4) return 0;
        g.add_edge(u, v, c, d);
        big += abs(d);
    }

    /// The return arc pays more than any path costs, so the cheapest circulation pushes as much as it can, capped at f
    int back = g.add_edge(n - 1, 0, f, -big);
    long long cost = g.solve();
    printf("%lld\n", g.flow(back) < f ? -1 : cost + big * f);
    return 0;
}
