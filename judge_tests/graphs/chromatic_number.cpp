// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/chromatic_number
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/chromatic_number.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    ChromaticNumber g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    printf("%d\n", g.solve());
    return 0;
}
