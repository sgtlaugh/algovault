// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/longest_increasing_subsequence
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/lis.cpp"
#undef main

static char buf[1 << 25];
int buf_len, buf_pos;

int read_int(){
    while (buf_pos < buf_len && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    int x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);

    int n = read_int();
    vector<int> a(n);
    for (auto& x : a) x = read_int();

    vector<int> res = lis_indices(a);

    string out = to_string(res.size()) + "\n";
    for (size_t i = 0; i < res.size(); i++){
        out += to_string(res[i]);
        out += i + 1 < res.size() ? ' ' : '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
