// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/inverse_matrix_mod_2
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<string> rows(n);
    for (auto& s : rows) cin >> s;

    vector<bitset<MAX>> inv;
    if (!matrix_inverse(n, rows_of(rows), inv)){
        cout << "-1\n";
        return 0;
    }

    string line(n, '0');
    for (const auto& row : inv){
        for (int j = 0; j < n; j++) line[j] = '0' + row[j];
        cout << line << '\n';
    }
    return 0;
}
