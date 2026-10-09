#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_2D.cpp"
#undef main

long long rect_sum(const vector<vector<long long>>& g, int i, int j, int k, int l){
    long long s = 0;
    for (int x = i; x <= k; x++) for (int y = j; y <= l; y++) s += g[x][y];
    return s;
}

/// All three variants against a plain grid, rectangles touching the borders included
int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 8 : 40), m = stress::rand_int(1, it % 10 ? 8 : 40);
        FenwickPointUpdate2D<long long> point(n, m);
        FenwickRangeUpdate2D<long long> range(n, m);
        FenwickFull2D<long long> full(n, m);
        vector<vector<long long>> a(n + 1, vector<long long>(m + 1)), b = a;

        for (int op = 0; op < 100; op++){
            int i = stress::rand_int(1, n), k = stress::rand_int(1, n), j = stress::rand_int(1, m), l = stress::rand_int(1, m);
            if (op % 5){  /// otherwise often empty, which must be a no-op
                if (i > k) swap(i, k);
                if (j > l) swap(j, l);
            }
            long long v = stress::rand_int(-1000000, 1000000);

            if (stress::rand_int(0, 1)){
                point.update(i, j, v);
                a[i][j] += v;
                range.update(i, j, k, l, v);
                full.update(i, j, k, l, v);
                for (int x = i; x <= k; x++) for (int y = j; y <= l; y++) b[x][y] += v;
            }
            else{
                assert(point.query(i, j, k, l) == rect_sum(a, i, j, k, l));
                assert(range.query(i, j) == b[i][j]);
                assert(full.query(i, j, k, l) == rect_sum(b, i, j, k, l));
            }
        }
    }

    return 0;
}
