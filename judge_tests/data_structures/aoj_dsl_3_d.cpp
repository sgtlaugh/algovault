// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_3_D
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/sliding_window_aggregation.cpp"
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
    int n = read_int(), k = read_int();

    MonotonicQueue<int> window;
    string out;
    for (int i = 0; i < n; i++){
        window.push(read_int());
        if (i >= k) window.pop();
        if (i >= k - 1){
            if (i > k - 1) out += ' ';
            out += to_string(window.best());
        }
    }

    out += '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
