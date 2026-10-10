// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/disjoint_sparse_table.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<long long> a(n);
    for (auto& x: a) x = read_int();

    auto table = DisjointST<long long>(a, 0);
    while (q--){
        int l = read_int(), r = read_int();
        printf("%lld\n", table.query(l, r - 1));
    }
    return 0;
}
