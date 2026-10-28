// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/chordal_graph_recognition
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/chordal_graph.cpp"
#undef main

static char buf[1 << 23];
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

    ChordalGraph g(n);
    for (int i = 0; i < m; i++){
        int u = read_int(), v = read_int();
        g.add_edge(u, v);
    }

    string out;
    if (g.run()){
        out = "YES\n";
        for (int v : g.order) out += to_string(v) + " ";
    }
    else{
        out = "NO\n" + to_string(g.cycle.size()) + "\n";
        for (int v : g.cycle) out += to_string(v) + " ";
    }
    out += "\n";
    fputs(out.c_str(), stdout);
    return 0;
}
