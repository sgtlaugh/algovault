// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_1_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/johnsons_algorithm.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    Johnson g(n);
    for (int i = 0; i < m; i++){
        int s, t;
        long long d;
        if (scanf("%d %d %lld", &s, &t, &d) != 3) return 0;
        g.add_edge(s, t, d);
    }

    if (!g.solve()){
        puts("NEGATIVE CYCLE");
        return 0;
    }

    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++){
            if (v) putchar(' ');
            if (g.dist[u][v] == Johnson::INF) printf("INF");
            else printf("%lld", g.dist[u][v]);
        }
        putchar('\n');
    }
    return 0;
}
