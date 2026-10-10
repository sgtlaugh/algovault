// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_add_rectangle_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_2D_sparse.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int rank_below(const vector<int>& v, int x){
    return lower_bound(v.begin(), v.end(), x) - v.begin();
}

int main(){
    int n = read_int(), q = read_int();
    vector<array<int, 5>> ops(n + q);
    for (int i = 0; i < n; i++) ops[i] = {0, read_int(), read_int(), read_int(), 0};
    for (int i = n; i < n + q; i++){
        if (read_int() == 0) ops[i] = {0, read_int(), read_int(), read_int(), 0};
        else ops[i] = {1, read_int(), read_int(), read_int(), read_int()};
    }

    vector<int> xs, ys;
    for (auto& op: ops){
        if (op[0] == 0) xs.push_back(op[1]), ys.push_back(op[2]);
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    /// 2e5 updates on 18 x 18 level chains create about 1.6e7 nodes, 2^25 slots keeps the load near half
    auto fen = FenwickSparse2D<long long>(xs.size(), ys.size(), 25);
    for (auto& op: ops){
        if (op[0] == 0) fen.update(rank_below(xs, op[1]) + 1, rank_below(ys, op[2]) + 1, op[3]);
        else{
            int l = rank_below(xs, op[1]) + 1, d = rank_below(ys, op[2]) + 1;
            int r = rank_below(xs, op[3]), u = rank_below(ys, op[4]);
            printf("%lld\n", fen.query(l, d, r, u));
        }
    }
    return 0;
}
