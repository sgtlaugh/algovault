// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/system_of_linear_equations_mod_2
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

/// 4096 variables plus the right-hand side bit need more than the default width
const size_t WIDTH = GF2_WIDTH + 64;

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

    auto [solvable, x, free_variables] = gauss(m, equations);
    if (!solvable){
        cout << "-1\n";
        return 0;
    }

    string line(m, '0');
    for (int j = 0; j < m; j++) line[j] = '0' + x[j];
    cout << free_variables << '\n' << line << '\n';

    /// gauss hides its reduced rows, eliminate_gf2 redoes them: free column f gives x_f = 1, x_j = bit f of row pos[j]
    vector<int> pos = eliminate_gf2(equations, m);
    for (int f = 0; f < m; f++){
        if (pos[f] != -1) continue;
        for (int j = 0; j < m; j++) line[j] = '0' + (j == f || (pos[j] != -1 && equations[pos[j]][f]));
        cout << line << '\n';
    }

    return 0;
}
