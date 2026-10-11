// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_count_distinct
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
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
    vector<int> values = a;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    for (int& x : a) x = lower_bound(values.begin(), values.end(), x) - values.begin();

    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        l = read_int(), r = read_int() - 1;  /// half-open [l, r) to inclusive, l = r gives an empty range
    }

    vector<int> cnt(values.size(), 0);
    int distinct = 0;
    auto add = [&](int i){ distinct += cnt[a[i]]++ == 0; };
    auto remove = [&](int i){ distinct -= --cnt[a[i]] == 0; };
    auto res = mo(n, queries, add, remove, [&](int){ return distinct; });

    string out;
    for (int x : res){
        out += to_string(x);
        out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
