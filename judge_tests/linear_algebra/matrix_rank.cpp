// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_rank
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_prime_mod.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<int>> a(n, vector<int>(m));
    for (auto& row : a){
        for (auto& x : row) scanf("%d", &x);
    }

    printf("%d\n", matrix_rank(a, 998244353));
    return 0;
}
