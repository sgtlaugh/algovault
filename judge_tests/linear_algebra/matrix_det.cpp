// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_det
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/determinant.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& x : row) scanf("%lld", &x);
    }

    printf("%lld\n", determinant_mod(a, 998244353));
    return 0;
}
