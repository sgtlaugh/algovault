// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/associative_array
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/hashmap.cpp"
#undef main

static char buf[1 << 26];
int buf_len = 0, buf_pos = 0;

long long read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    long long x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int q = read_int();

    auto hashmap = HashMap<long long, long long>(q);
    string out;
    while (q--){
        int t = read_int();
        long long k = read_int();
        if (t == 0) hashmap.set(k, read_int());
        else{
            out += to_string(hashmap.get(k));
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
