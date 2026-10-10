// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_mode_query
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<int> a(n);
    for (auto& x: a) x = read_int();

    vector<int> values = a;
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    vector<int> code(n);
    for (int i = 0; i < n; i++) code[i] = lower_bound(values.begin(), values.end(), a[i]) - values.begin();

    vector<pair<int, int>> queries(q);
    for (auto& [l, r]: queries){
        l = read_int(), r = read_int() - 1;
    }

    /// (count, code) of the current mode, log keeps every incremented code so rollback can undo them
    vector<int> cnt(values.size(), 0), log;
    pair<int, int> best = {0, 0};
    vector<pair<int, pair<int, int>>> saved;
    vector<pair<int, int>> res(q);

    auto add = [&](int i){
        int c = code[i];
        log.push_back(c);
        best = max(best, {++cnt[c], c});
    };
    auto snapshot = [&](){
        saved.push_back({(int)log.size(), best});
    };
    auto rollback = [&](){
        auto [size, prev_best] = saved.back();
        saved.pop_back();
        for (; (int)log.size() > size; log.pop_back()) cnt[log.back()]--;
        best = prev_best;
    };
    auto answer = [&](int qi){
        res[qi] = best;
    };
    mo_with_rollback(n, queries, add, snapshot, rollback, answer);

    for (auto& [count, c]: res) printf("%d %d\n", values[c], count);
    return 0;
}
