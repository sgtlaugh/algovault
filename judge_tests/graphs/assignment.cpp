// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/assignment
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/hungarian_algorithm.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& x : row){
            if (scanf("%lld", &x) != 1) return 0;
        }
    }

    auto [cost, pairs] = hungarian<long long>(a, true);
    vector<int> p(n);
    for (auto [row, col] : pairs) p[row] = col;

    printf("%lld\n", cost);
    for (int i = 0; i < n; i++) printf("%d%c", p[i], i + 1 < n ? ' ' : '\n');
    return 0;
}
