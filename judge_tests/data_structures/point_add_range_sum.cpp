// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_range_sum
// competitive-verifier: TLE 5
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

    vector<long long> a(n);
    for (auto& x : a) x = read_int();
    auto fen = FenwickPointUpdate<long long>(a);

    string out;
    while (q--){
        int t = read_int();
        if (t == 0){
            int p = read_int();
            fen.update(p, read_int());
        }
        else{
            int l = read_int(), r = read_int();
            out += to_string(fen.query(l, r - 1));
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
