// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/predecessor_problem
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fast_set.cpp"
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

    FastSet fs(n);
    for (int i = 0; i < n; i++, c = read_char()){
        if (c == '1') fs.insert(i);
    }

    string out;
    while (q--){
        int type = read_int(), k = read_int();
        if (type == 0) fs.insert(k);
        if (type == 1) fs.erase(k);
        if (type == 2) out += fs.contains(k) ? "1\n" : "0\n";
        if (type == 3){
            int x = fs.next(k);
            out += to_string(x < n ? x : -1) + '\n';
        }
        if (type == 4) out += to_string(fs.prev(k)) + '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
