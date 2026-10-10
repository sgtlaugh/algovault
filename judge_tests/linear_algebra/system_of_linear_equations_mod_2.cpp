// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/system_of_linear_equations_mod_2
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

/// 4096 variables plus the right-hand side bit need more than the default width
const size_t WIDTH = MAX + 64;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if (!(cin >> n >> m)) return 0;
    vector<bitset<WIDTH>> equations(n);
    string s;
    for (int i = 0; i < n; i++){
        cin >> s;
        for (int j = 0; j < m; j++) equations[i][j] = s[j] == '1';
    }
    cin >> s;
    for (int i = 0; i < n; i++) equations[i][m] = s[i] == '1';

    bitset<WIDTH> res;
    int f_var = gauss(m, equations, res);
    if (f_var == -1){
        cout << "-1\n";
        return 0;
    }

    string line(m, '0');
    for (int j = 0; j < m; j++) line[j] = '0' + res[j];
    cout << f_var << '\n' << line << '\n';

    /// gauss does not return its reduced rows, so they are rebuilt with eliminate_gf2
    /// Each free column f gives a kernel vector: x_f = 1, and x_j = bit f of the pivot row of each pivot column j
    vector<int> pos = eliminate_gf2(equations, m);
    for (int f = 0; f < m; f++){
        if (pos[f] != -1) continue;
        for (int j = 0; j < m; j++) line[j] = '0' + (j == f || (pos[j] != -1 && equations[pos[j]][f]));
        cout << line << '\n';
    }
    return 0;
}
