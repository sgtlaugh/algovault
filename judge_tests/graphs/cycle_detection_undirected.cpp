// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/cycle_detection_undirected
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/find_cycle.cpp"
#undef main

static char buf[1 << 24];
int buf_len = 0, buf_pos = 0;

inline int read_int(){
    while (buf[buf_pos] < '0') buf_pos++;
    int x = 0;
    while (buf[buf_pos] >= '0') x = x * 10 + buf[buf_pos++] - '0';
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf) - 1, stdin);
    buf[buf_len] = 0;
    int n = read_int(), m = read_int();

    FindCycle g(n, FindCycle::UNDIRECTED);
    for (int i = 0; i < m; i++){
        int u = read_int(), v = read_int();
        g.add_edge(u, v);
    }

    FindCycle::Cycle c = g.solve();
    if (!c.exists){
        puts("-1");
        return 0;
    }

    string out = to_string(c.vertices.size()) + "\n";
    for (size_t i = 0; i < c.vertices.size(); i++) out += to_string(c.vertices[i]) + (i + 1 < c.vertices.size() ? " " : "\n");
    for (size_t i = 0; i < c.edges.size(); i++) out += to_string(c.edges[i]) + (i + 1 < c.edges.size() ? " " : "\n");
    fputs(out.c_str(), stdout);
    return 0;
}
