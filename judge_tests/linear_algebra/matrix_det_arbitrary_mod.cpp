// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_det_arbitrary_mod
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/determinant.cpp"
#undef main

int main(){
    int n;
    long long m;
    if (scanf("%d %lld", &n, &m) != 2) return 0;
    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& x : row) scanf("%lld", &x);
    }

    printf("%lld\n", determinant_mod(a, m));
    return 0;
}
