// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_kth_smallest
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/persistent_segment_tree.cpp"
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

    vector<int> a(n);
    for (int& x : a) x = read_int();
    auto rk = RangeKth<int>(a);

    string out;
    while (q--){
        int l = read_int(), r = read_int(), k = read_int();
        out += to_string(rk.kth(l, r - 1, k));
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
