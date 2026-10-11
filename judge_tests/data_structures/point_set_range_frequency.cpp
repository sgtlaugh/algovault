// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_frequency
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/mo_with_updates.cpp"
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

    vector<pair<int, int>> updates;
    vector<array<int, 3>> queries;
    vector<int> targets;
    for (int i = 0; i < q; i++){
        if (read_int() == 0){
            int p = read_int(), x = read_int();
            updates.push_back({p, x});
            continue;
        }

        int l = read_int(), r = read_int(), x = read_int();
        queries.push_back({l, r - 1, (int)updates.size()}), targets.push_back(x);  /// half-open [l, r) to inclusive
    }

    vector<int> values = a;
    for (auto& [p, x] : updates) values.push_back(x);
    values.insert(values.end(), targets.begin(), targets.end());
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    auto compress = [&](int& x){ x = lower_bound(values.begin(), values.end(), x) - values.begin(); };
    for (int& x : a) compress(x);
    for (auto& [p, x] : updates) compress(x);
    for (int& x : targets) compress(x);

    vector<int> cnt(values.size(), 0);
    auto add = [&](int i){ cnt[a[i]]++; };
    auto remove = [&](int i){ cnt[a[i]]--; };
    auto res = mo_with_updates(a, queries, updates, add, remove, [&](int qi){ return cnt[targets[qi]]; });

    string out;
    for (int x : res){
        out += to_string(x);
        out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
