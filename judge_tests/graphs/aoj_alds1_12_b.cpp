// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_12_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/dijkstra.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    DenseDijkstra g(n);
    for (int i = 0; i < n; i++){
        int u, k;
        if (scanf("%d %d", &u, &k) != 2) return 0;
        while (k--){
            int v;
            long long c;
            if (scanf("%d %lld", &v, &c) != 2) return 0;
            g.add_edge(u, v, c);
        }
    }

    g.run(0);
    for (int v = 0; v < n; v++) printf("%d %lld\n", v, g.dist[v]);
    return 0;
}
