#include "common.h"

#define main library_main
#include "../code_library/fenwick_tree_2D_implicit.cpp"
#undef main

/// Rectangle updates and point queries
int main(){
    /// Small grids against a full brute force array, rectangles often touch row and column n
    for (long long round = 0; round < stress::scaled(20); round++){
        const int n = stress::rand_int(1, 30);
        FenwickImplicit2D<long long> small(n);
        vector<vector<long long>> grid(n + 1, vector<long long>(n + 1, 0));
        for (int it = 0; it < 150; it++){
            int x1 = stress::rand_int(1, n), x2 = stress::rand_int(1, n), y1 = stress::rand_int(1, n), y2 = stress::rand_int(1, n);
            if (x1 > x2) swap(x1, x2);
            if (y1 > y2) swap(y1, y2);
            if (stress::rand_int(0, 3) == 0) x2 = y2 = n;
            long long v = stress::rand_int(-100, 100);
            small.update(x1, y1, x2, y2, v);
            for (int x = x1; x <= x2; x++) for (int y = y1; y <= y2; y++) grid[x][y] += v;

            int x = stress::rand_int(1, n), y = stress::rand_int(1, n);
            assert(small.query(x, y) == grid[x][y]);
        }
        for (int x = 1; x <= n; x++) for (int y = 1; y <= n; y++) assert(small.query(x, y) == grid[x][y]);
    }

    /// Large grid against the list of updates, at most 1500 updates (~0.8M pool nodes) at any scale, larger scales add queries
    const int N = 1000000;
    FenwickImplicit2D<long long> large(N);
    vector<array<long long, 5>> rects;
    for (long long it = 0; it < stress::scaled(1500); it++){
        int x1 = stress::rand_int(1, N), x2 = stress::rand_int(1, N), y1 = stress::rand_int(1, N), y2 = stress::rand_int(1, N);
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);
        if (rects.size() < 1500){
            long long v = stress::rand_int(-100, 100);
            large.update(x1, y1, x2, y2, v);
            rects.push_back({x1, y1, x2, y2, v});
        }

        int x = stress::rand_int(0, 1) ? stress::rand_int(1, N) : x1, y = stress::rand_int(0, 1) ? stress::rand_int(1, N) : y2;
        long long expected = 0;
        for (auto& r : rects) if (r[0] <= x && x <= r[2] && r[1] <= y && y <= r[3]) expected += r[4];
        assert(large.query(x, y) == expected);
    }
    return 0;
}
