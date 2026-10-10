// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/minimum_spanning_tree
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/minimum_spanning_tree.cpp"
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
    int n = read_int(), m = read_int();
    Kruskal kruskal(n);
    Boruvka boruvka(n);
    for (int i = 0; i < m; i++){
        int a = read_int(), b = read_int();
        long long c = read_int();
        kruskal.add_edge(a, b, c);
        boruvka.add_edge(a, b, c);
    }

    auto res = kruskal.solve();
    assert(boruvka.solve().weight == res.weight);

    string out = to_string(res.weight) + '\n';
    for (int i = 0; i < (int)res.chosen.size(); i++){
        if (i) out += ' ';
        out += to_string(res.chosen[i]);
    }
    out += '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
