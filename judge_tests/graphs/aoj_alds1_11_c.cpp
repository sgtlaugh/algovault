// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_11_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/dijkstra.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    ZeroOneBFS g(n);
    for (int i = 0; i < n; i++){
        int u, k;
        if (scanf("%d %d", &u, &k) != 2) return 0;
        for (int j = 0; j < k; j++){
            int v;
            if (scanf("%d", &v) != 1) return 0;
            g.add_edge(u - 1, v - 1, 1);
        }
    }

    g.run(0);
    for (int v = 0; v < n; v++) printf("%d %lld\n", v + 1, g.dist[v] == DIJKSTRA_INF ? -1 : g.dist[v]);
    return 0;
}
