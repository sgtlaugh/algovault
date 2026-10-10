// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_frequency
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

    vector<int> update_pos, value, targets, res(q, -1), index;
    vector<array<int, 3>> queries;
    for (int i = 0; i < q; i++){
        if (read_int() == 0){
            update_pos.push_back(read_int());
            value.push_back(read_int());
            continue;
        }

        int l = read_int(), r = read_int(), x = read_int();
        res[i] = 0;
        if (l < r) queries.push_back({l, r - 1, (int)update_pos.size()}), targets.push_back(x), index.push_back(i);
    }

    vector<int> values = a;
    values.insert(values.end(), value.begin(), value.end());
    values.insert(values.end(), targets.begin(), targets.end());
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    auto compress = [&](int& x){ x = lower_bound(values.begin(), values.end(), x) - values.begin(); };
    for (int& x : a) compress(x);
    for (int& x : value) compress(x);
    for (int& x : targets) compress(x);

    vector<int> cnt(values.size(), 0);
    auto add = [&](int i){ cnt[a[i]]++; };
    auto remove = [&](int i){ cnt[a[i]]--; };
    auto toggle = [&](int j){ swap(a[update_pos[j]], value[j]); };
    mo_with_updates(n, queries, update_pos, add, remove, toggle, [&](int qi){ res[index[qi]] = cnt[targets[qi]]; });

    string out;
    for (int x : res){
        if (x == -1) continue;
        out += to_string(x);
        out += '\n';
    }
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
