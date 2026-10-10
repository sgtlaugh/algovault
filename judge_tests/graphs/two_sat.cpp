// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/two_sat
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/2SAT_kosaraju.cpp"
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

/// Skips everything that cannot start a number, which also skips the "p cnf" prefix
int read_int(){
    int c = read_char();
    while (c != -1 && c != '-' && (c < '0' || c > '9')) c = read_char();
    bool neg = c == '-';
    if (neg) c = read_char();

    int x = 0;
    while (c >= '0' && c <= '9') x = x * 10 + (c - '0'), c = read_char();
    return neg ? -x : x;
}

int main(){
    int n = read_int(), m = read_int();
    Graph g(n);
    for (int i = 0; i < m; i++){
        int a = read_int(), b = read_int();
        read_int();
        g.add_or(a, b);
    }

    if (!g.is_satisfiable()){
        puts("s UNSATISFIABLE");
        return 0;
    }

    string out = "s SATISFIABLE\nv";
    for (int i = 1; i <= n; i++) out += ' ' + to_string(g.value(i) ? i : -i);
    out += " 0\n";
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
