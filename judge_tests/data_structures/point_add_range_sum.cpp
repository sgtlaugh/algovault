// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

long long read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    long long x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    auto fen = FenwickPointUpdate<long long>(n);
    for (int i = 1; i <= n; i++) fen.update(i, read_int());

    string out;
    while (q--){
        int t = read_int();
        if (t == 0){
            int p = read_int();
            fen.update(p + 1, read_int());
        }
        else{
            int l = read_int(), r = read_int();
            out += to_string(fen.query(l + 1, r));
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
