// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/inverse_matrix
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_prime_mod.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<vector<int>> a(n, vector<int>(n));
    for (auto& row : a){
        for (auto& x : row) scanf("%d", &x);
    }

    vector<vector<int>> inv;
    if (!matrix_inverse(a, inv, 998244353)){
        puts("-1");
        return 0;
    }
    for (auto& row : inv){
        for (int j = 0; j < n; j++) printf("%d%c", row[j], j + 1 == n ? '\n' : ' ');
    }
    return 0;
}
