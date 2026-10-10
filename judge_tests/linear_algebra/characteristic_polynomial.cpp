// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/characteristic_polynomial
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/characteristic_polynomial.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& x : row) scanf("%lld", &x);
    }

    vector<long long> c = characteristic_polynomial(a, 998244353);

    for (int i = 0; i <= n; i++) printf("%lld%c", c[i], i == n ? '\n' : ' ');
    return 0;
}
