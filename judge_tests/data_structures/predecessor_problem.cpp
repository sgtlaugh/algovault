// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/predecessor_problem
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

static char buffer[1 << 16];
static size_t buffer_len = 0, buffer_pos = 0;

int read_char(){
    if (buffer_pos == buffer_len){
        buffer_len = fread(buffer, 1, sizeof(buffer), stdin);
        buffer_pos = 0;
        if (buffer_len == 0) return EOF;
    }
    return buffer[buffer_pos++];
}

int read_int(){
    int c = read_char();
    while (c == ' ' || c == '\n' || c == '\r') c = read_char();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = read_char()) x = x * 10 + (c - '0');
    return x;
}

int main(){
    int n = read_int(), q = read_int();

    int c = read_char();
    while (c != '0' && c != '1') c = read_char();

    /// Key k lives at Fenwick index k + 1
    vector<char> present(n, 0);
    FenwickPointUpdate<int> fenwick(n);
    int total = 0;
    for (int i = 0; i < n; i++, c = read_char()){
        if (c == '1') present[i] = 1, fenwick.update(i + 1, 1), total++;
    }

    string out;
    while (q--){
        int type = read_int(), k = read_int();
        if (type == 0 && !present[k]) present[k] = 1, fenwick.update(k + 1, 1), total++;
        if (type == 1 && present[k]) present[k] = 0, fenwick.update(k + 1, -1), total--;
        if (type == 2) out += present[k] ? "1\n" : "0\n";
        if (type == 3){
            int before = fenwick.query(k);
            out += to_string(before < total ? fenwick.lower_bound(before + 1) - 1 : -1) + '\n';
        }
        if (type == 4){
            int upto = fenwick.query(k + 1);
            out += to_string(upto > 0 ? fenwick.lower_bound(upto) - 1 : -1) + '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
