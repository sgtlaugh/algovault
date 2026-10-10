// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/double_ended_priority_queue
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/treap.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

int read_int(){
    while (buf_pos < buf_len && buf[buf_pos] != '-' && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    bool negative = buf_pos < buf_len && buf[buf_pos] == '-';
    if (negative) buf_pos++;

    int x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return negative ? -x : x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    Treap<int> treap;
    for (int i = 0; i < n; i++) treap.insert(read_int());

    string out;
    while (q--){
        int t = read_int();
        if (t == 0){
            treap.insert(read_int());
            continue;
        }

        int x = treap.kth(t == 1 ? 0 : treap.size() - 1);
        treap.erase(x);
        out += to_string(x);
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
