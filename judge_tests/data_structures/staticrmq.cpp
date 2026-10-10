// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/staticrmq
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/sparse_table.cpp"
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

    auto table = SparseTable<int>(a);
    LinearSparseTable<int> linear(a);
    while (q--){
        int l = read_int(), r = read_int() - 1;
        int res = table.query(l, r);
        assert(res == linear.query(l, r));
        printf("%d\n", res);
    }
    return 0;
}
