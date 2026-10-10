// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/pow_of_matrix
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/matrix.cpp"
#undef main

int main(){
    int n;
    long long k;
    if (scanf("%d %lld", &n, &k) != 2) return 0;
    Matrix a(n, n, 998244353);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) scanf("%lld", &a[i][j]);
    }

    Matrix b = a.pow(k);

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) printf("%lld%c", b[i][j], j + 1 == n ? '\n' : ' ');
    }
    return 0;
}
