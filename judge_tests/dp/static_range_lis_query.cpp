// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_lis_query
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/range_lis.cpp"
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

    int n = read_int(), q = read_int();
    vector<int> p(n);
    for (auto& x : p) x = read_int();

    /// The problem's ranges are half-open [l, r)
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries) l = read_int(), r = read_int() - 1;

    RangeLIS lis(p);
    vector<int> res = lis.query(queries);

    string out;
    for (int x : res) out += to_string(x) + '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
