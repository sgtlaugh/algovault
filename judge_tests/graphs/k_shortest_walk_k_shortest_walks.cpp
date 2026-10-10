// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/k_shortest_walk
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/k_shortest_walks.cpp"
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
    int n = read_int(), m = read_int(), s = read_int(), t = read_int(), k = read_int();

    KShortestWalks g(n);
    for (int i = 0; i < m; i++){
        int u = read_int(), v = read_int(), c = read_int();
        g.add_edge(u, v, c);
    }

    auto lengths = g.solve(s, t, k);
    lengths.resize(k, -1);

    string out;
    for (long long x : lengths) out += to_string(x) + "\n";
    fputs(out.c_str(), stdout);
    return 0;
}
