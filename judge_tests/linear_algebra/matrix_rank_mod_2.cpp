// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/matrix_rank_mod_2
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_bitset.cpp"
#undef main

/// Rows go along the longer side, width = smallest power of two >= min(n, m): 2^24 rows of 4096 bits would take 8 GB
template <size_t W>
int read_rank(int n, int m){
    if constexpr (W < GF2_WIDTH){
        if ((int)W < min(n, m)) return read_rank<2 * W>(n, m);
    }

    bool transpose = m > n;
    vector<bitset<W>> a(max(n, m));
    string s;
    for (int i = 0; i < n; i++){
        cin >> s;
        for (int j = 0; j < m; j++){
            if (s[j] == '0') continue;
            if (transpose) a[j][i] = 1;
            else a[i][j] = 1;
        }
    }

    return matrix_rank(min(n, m), move(a));
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    if (!(cin >> n >> m)) return 0;
    cout << read_rank<64>(n, m) << '\n';
    return 0;
}
