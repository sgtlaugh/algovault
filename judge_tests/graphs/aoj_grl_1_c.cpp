// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_1_C
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/floyd_warshall.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    FloydWarshall fw(n);
    for (int i = 0; i < m; i++){
        int s, t;
        long long d;
        if (scanf("%d %d %lld", &s, &t, &d) != 3) return 0;
        fw.add_edge(s, t, d);
    }

    fw.solve();
    for (int v = 0; v < n; v++){
        if (fw.dist[v][v] == FloydWarshall::NEG_INF){
            puts("NEGATIVE CYCLE");
            return 0;
        }
    }

    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++){
            if (v) putchar(' ');
            if (fw.dist[u][v] == FloydWarshall::INF) printf("INF");
            else printf("%lld", fw.dist[u][v]);
        }
        putchar('\n');
    }
    return 0;
}
