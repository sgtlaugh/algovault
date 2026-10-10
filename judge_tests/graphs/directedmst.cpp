// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/directedmst
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/directed_mst.cpp"
#undef main

int main(){
    int n, m, s;
    if (scanf("%d %d %d", &n, &m, &s) != 3) return 0;

    DirectedMST g(n);
    for (int i = 0; i < m; i++){
        int u, v;
        long long w;
        if (scanf("%d %d %lld", &u, &v, &w) != 3) return 0;
        g.add_edge(u, v, w);
    }

    auto res = g.solve(s);
    res.parent[s] = s;

    string out = to_string(res.weight) + "\n";
    for (int v = 0; v < n; v++) out += to_string(res.parent[v]) + (v + 1 < n ? " " : "\n");
    fputs(out.c_str(), stdout);
    return 0;
}
