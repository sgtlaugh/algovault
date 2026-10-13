#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fenwick_tree_2D_implicit.cpp"
#undef main

/// Rectangle updates and point queries
int main(){
    /// Small grids against a full brute force array, square or not, rectangles often touch the last row and column
    /// Few rounds since every instance reserves the whole node pool
    for (long long round = 0; round < stress::scaled(40); round++){
        const int n = stress::rand_int(1, 30), m = round % 2 ? n : stress::rand_int(1, 30);
        FenwickImplicit2D<long long> small(n, m);
        vector<vector<long long>> grid(n, vector<long long>(m, 0));
        for (int it = 0; it < 150; it++){
            int x1 = stress::rand_int(0, n - 1), x2 = stress::rand_int(0, n - 1), y1 = stress::rand_int(0, m - 1), y2 = stress::rand_int(0, m - 1);
            if (x1 > x2) swap(x1, x2);
            if (y1 > y2) swap(y1, y2);
            if (stress::rand_int(0, 3) == 0) x2 = n - 1, y2 = m - 1;
            if (it % 7 == 0) swap(x1, x2);  /// an empty rectangle when x1 != x2, the update must change nothing

            long long v = stress::rand_int(-100, 100);
            small.update(x1, y1, x2, y2, v);
            for (int x = x1; x <= x2; x++) for (int y = y1; y <= y2; y++) grid[x][y] += v;

            int x = stress::rand_int(0, n - 1), y = stress::rand_int(0, m - 1);
            assert(small.query(x, y) == grid[x][y]);
        }

        for (int x = 0; x < n; x++) for (int y = 0; y < m; y++) assert(small.query(x, y) == grid[x][y]);
    }

    /// Large grid against the list of updates, at most 1500 updates (~1.2M pool nodes) at any scale, larger scales add queries
    /// 2e9 columns, so a column midpoint written as (a + b) / 2 overflows int
    const int N = 1000000, M = 2000000000;
    FenwickImplicit2D<long long> large(N, M);
    vector<array<long long, 5>> rects;
    for (long long it = 0; it < stress::scaled(1500); it++){
        int x1 = stress::rand_int(0, N - 1), x2 = stress::rand_int(0, N - 1), y1 = stress::rand_int(0, M - 1), y2 = stress::rand_int(0, M - 1);
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);

        if (rects.size() < 1500){
            long long v = stress::rand_int(-100, 100);
            large.update(x1, y1, x2, y2, v);
            rects.push_back({x1, y1, x2, y2, v});
        }

        int x = stress::rand_int(0, 1) ? stress::rand_int(0, N - 1) : x1, y = stress::rand_int(0, 1) ? stress::rand_int(0, M - 1) : y2;
        long long expected = 0;
        for (auto& r : rects) if (r[0] <= x && x <= r[2] && r[1] <= y && y <= r[3]) expected += r[4];
        assert(large.query(x, y) == expected);
    }

    return 0;
}
