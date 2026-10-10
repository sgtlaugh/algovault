// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/bipartitematching
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/maxflow.cpp"
#undef main

int main(){
    int l, r, m;
    if (scanf("%d %d %d", &l, &r, &m) != 3) return 0;

    int src = l + r, sink = l + r + 1;
    FlowGraph g(l + r + 2, src, sink);
    for (int u = 0; u < l; u++) g.add_directed_edge(src, u, 1);
    for (int v = 0; v < r; v++) g.add_directed_edge(l + v, sink, 1);

    vector<int> ids(m);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        ids[i] = g.E.size();
        g.add_directed_edge(u, l + v, 1);
    }

    long long k = g.maxflow();

    string out = to_string(k) + "\n";
    for (int id : ids){
        const Edge& e = g.E[id];
        if (e.flow == 1) out += to_string(e.u) + " " + to_string(e.v - l) + "\n";
    }
    fputs(out.c_str(), stdout);
    return 0;
}
