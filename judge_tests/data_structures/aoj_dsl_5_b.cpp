// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_5_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_2D.cpp"
#undef main

int main(){
    const int size = 1000;
    int n;
    if (scanf("%d", &n) != 1) return 0;

    /// Unit cell [x, x + 1] x [y, y + 1] is Fenwick cell (x, y)
    FenwickRangeUpdate2D<int> fenwick(size, size);
    for (int i = 0; i < n; i++){
        int x1, y1, x2, y2;
        if (scanf("%d %d %d %d", &x1, &y1, &x2, &y2) != 4) return 0;
        fenwick.update(x1, y1, x2 - 1, y2 - 1, 1);
    }

    int best = 0;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++) best = max(best, fenwick.query(i, j));
    }

    printf("%d\n", best);
    return 0;
}
