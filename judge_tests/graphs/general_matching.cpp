// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/general_matching
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/graph_matching.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    Blossom b(n);
    Graph tutte(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        b.add_edge(u, v);
        tutte.add_edge(u, v);
    }

    int size = b.maximum_matching();
    assert(tutte.maximum_matching() == size);

    printf("%d\n", size);
    for (int v = 0; v < n; v++){
        if (b.mate[v] > v) printf("%d %d\n", v, b.mate[v]);
    }
    return 0;
}
