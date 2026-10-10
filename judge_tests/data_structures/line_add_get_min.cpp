// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/line_add_get_min
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/li_chao_tree.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

long long read_int(){
    while (buf_pos < buf_len && buf[buf_pos] != '-' && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    bool negative = buf_pos < buf_len && buf[buf_pos] == '-';
    if (negative) buf_pos++;

    long long x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return negative ? -x : x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    const long long C = 1000000000;
    LiChaoTree<long long> tree(-C, C);
    PersistentLiChaoTree<long long> persistent(-C, C);
    int version = 0;
    for (int i = 0; i < n; i++){
        long long a = read_int(), b = read_int();
        tree.add_line(a, b);
        version = persistent.add_line(version, a, b);
    }

    string out;
    while (q--){
        int t = read_int();
        if (t == 0){
            long long a = read_int(), b = read_int();
            tree.add_line(a, b);
            version = persistent.add_line(version, a, b);
        }
        else{
            long long p = read_int(), res = tree.query(p);
            assert(persistent.query(version, p) == res);
            out += to_string(res);
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
