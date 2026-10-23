// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_13_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/dancing_links.cpp"
#undef main

int main(){
    int k;
    if (scanf("%d", &k) != 1) return 0;

    vector<int> fixed_column(8, -1);
    for (int i = 0; i < k; i++){
        int r, c;
        if (scanf("%d %d", &r, &c) != 2) return 0;
        fixed_column[r] = c;
    }

    /// Columns 0-7 rows, 8-15 columns, 16-30 diagonals, 31-45 anti-diagonals, each diagonal with a slack row
    DancingLinks dlx(46);
    for (int r = 0; r < 8; r++){
        for (int c = 0; c < 8; c++){
            if (fixed_column[r] != -1 && fixed_column[r] != c) continue;
            dlx.add_row(r * 8 + c, {r, c + 8, r + c + 16, r - c + 38});
        }
    }
    for (int d = 16; d <= 45; d++) dlx.add_row(64 + d, {d});

    auto rows = dlx.exact_cover();
    if (!rows) return 0;

    vector<string> board(8, string(8, '.'));
    for (int id : *rows){
        if (id < 64) board[id / 8][id % 8] = 'Q';
    }
    for (auto& line : board) puts(line.c_str());
    return 0;
}
