// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/tree_path_composite_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/rerooting.cpp"
#undef main

const long long MOD = 998244353;

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int();
    vector<long long> a(n);
    for (auto& x : a) x = read_int();
    using Value = pair<long long, long long>;
    Rerooting<Value, pair<long long, long long>> tree(n);
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int(), b = read_int(), c = read_int();
        tree.add_edge(u, v, {b, c});
    }

    /// A value is (nodes strictly below, sum of P over them), so from's own a joins as the value crosses the edge
    auto merge = [](const Value& x, const Value& y){
        return Value{x.first + y.first, (x.second + y.second) % MOD};
    };
    auto apply_edge = [&](const Value& x, int from, int, const pair<long long, long long>& e){
        long long cnt = x.first + 1;
        return Value{cnt, (e.first * ((x.second + a[from]) % MOD) + e.second * (cnt % MOD)) % MOD};
    };
    vector<Value> res = tree.solve(merge, Value{0, 0}, apply_edge);

    string out;
    for (int x = 0; x < n; x++) out += to_string((res[x].second + a[x]) % MOD) + (x + 1 < n ? ' ' : '\n');
    fputs(out.c_str(), stdout);
    return 0;
}
