// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lca
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/link_cut_tree.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

int read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    int x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    auto lct = LinkCutTree<int>(n);
    for (int i = 1; i < n; i++) lct.link(i, read_int());

    string out;
    while (q--){
        int u = read_int(), v = read_int();
        out += to_string(lct.lca(u, v));
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
