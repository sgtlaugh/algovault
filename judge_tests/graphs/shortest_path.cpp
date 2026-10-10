// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/shortest_path
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/dijkstra.cpp"
#undef main

static char in_buf[1 << 16];
static size_t in_len = 0, in_pos = 0;

int read_char(){
    if (in_pos == in_len){
        in_len = fread(in_buf, 1, sizeof(in_buf), stdin), in_pos = 0;
        if (in_len == 0) return -1;
    }
    return in_buf[in_pos++];
}

long long read_int(){
    int c = read_char();
    while (c != -1 && (c < '0' || c > '9')) c = read_char();

    long long x = 0;
    while (c >= '0' && c <= '9') x = x * 10 + (c - '0'), c = read_char();
    return x;
}

int main(){
    int n = read_int(), m = read_int(), s = read_int(), t = read_int();
    Dijkstra g(n);
    for (int i = 0; i < m; i++){
        int a = read_int(), b = read_int();
        long long c = read_int();
        g.add_edge(a, b, c);
    }

    g.run(s);
    if (g.dist[t] == DIJKSTRA_INF){
        puts("-1");
        return 0;
    }

    auto path = g.path(t);
    string out = to_string(g.dist[t]) + ' ' + to_string(path.size() - 1) + '\n';
    for (int i = 0; i + 1 < (int)path.size(); i++) out += to_string(path[i]) + ' ' + to_string(path[i + 1]) + '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
