// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/maximum_independent_set
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/max_clique.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;

    CliqueGraph g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        g.add_edge(u, v);
    }

    auto set = g.max_independent_set();
    int x = set.size();

    printf("%d\n", x);
    for (int i = 0; i < x; i++) printf("%d%c", set[i], i + 1 < x ? ' ' : '\n');
    return 0;
}
