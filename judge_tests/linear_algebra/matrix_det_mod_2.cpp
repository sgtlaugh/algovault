// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_det_mod_2
// competitive-verifier: TLE 5
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

    /// The determinant is 1 mod 2 exactly when the matrix has full rank
    cout << (matrix_rank(n, rows_of(rows)) == n) << '\n';
    return 0;
}
