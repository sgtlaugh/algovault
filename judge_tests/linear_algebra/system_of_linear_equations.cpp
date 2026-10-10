// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/system_of_linear_equations
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_prime_mod.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<int>> equations(n, vector<int>(m + 1));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) scanf("%d", &equations[i][j]);
    }
    for (int i = 0; i < n; i++) scanf("%d", &equations[i][m]);

    vector<int> res;
    vector<vector<int>> basis;
    int r = gauss(equations, res, basis, 998244353);

    if (r == -1){
        puts("-1");
        return 0;
    }
    printf("%d\n", r);
    for (int j = 0; j < m; j++) printf("%d%c", res[j], j + 1 == m ? '\n' : ' ');
    for (auto& x : basis){
        for (int j = 0; j < m; j++) printf("%d%c", x[j], j + 1 == m ? '\n' : ' ');
    }
    return 0;
}
