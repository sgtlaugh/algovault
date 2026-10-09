#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sparse_table_2D.cpp"
#undef main

template <typename T>
vector<vector<T>> random_grid(int n, int m, T lo, T hi){
    vector<vector<T>> grid(n, vector<T>(m));
    for (auto& row : grid) for (auto& x : row) x = stress::rand_int(lo, hi);
    return grid;
}

template <typename T>
void check(const vector<vector<T>>& grid, const SparseTable2D<T>& mn, const SparseTable2D<T, greater<T>>& mx, int x1, int y1, int x2, int y2){
    T lo = grid[x1][y1], hi = grid[x1][y1];
    for (int i = x1; i <= x2; i++){
        for (int j = y1; j <= y2; j++) lo = min(lo, grid[i][j]), hi = max(hi, grid[i][j]);
    }
    assert(mn.query(x1, y1, x2, y2) == lo && mx.query(x1, y1, x2, y2) == hi);
}

/// Rectangle min and max against a direct scan, every rectangle for small grids and random ones for large grids
int main(){
    for (long long it = 0; it < stress::scaled(600); it++){
        int n = it < 144 ? it / 12 + 1 : stress::rand_int(1, 120);
        int m = it < 144 ? it % 12 + 1 : stress::rand_int(1, 120);
        auto grid = it % 3 == 0 ? random_grid<int>(n, m, -2, 2) : random_grid<int>(n, m, INT_MIN, INT_MAX);
        SparseTable2D<int> mn(grid);
        SparseTable2D<int, greater<int>> mx(grid);

        if (n * m <= 144){
            for (int x1 = 0; x1 < n; x1++) for (int x2 = x1; x2 < n; x2++){
                for (int y1 = 0; y1 < m; y1++) for (int y2 = y1; y2 < m; y2++) check(grid, mn, mx, x1, y1, x2, y2);
            }
            continue;
        }

        for (int q = 0; q < 300; q++){
            int x1 = stress::rand_int(0, n - 1), x2 = stress::rand_int(x1, n - 1);
            int y1 = stress::rand_int(0, m - 1), y2 = stress::rand_int(y1, m - 1);
            check(grid, mn, mx, x1, y1, x2, y2);
        }
    }

    /// long long extremes on a non-square grid with power of two and off-by-one sides
    for (auto [n, m] : vector<pair<int, int>>{{64, 65}, {127, 128}, {1, 300}, {300, 1}}){
        auto grid = random_grid<long long>(n, m, LLONG_MIN, LLONG_MAX);
        SparseTable2D<long long> mn(grid);
        SparseTable2D<long long, greater<long long>> mx(grid);
        for (int q = 0; q < 500; q++){
            int x1 = stress::rand_int(0, n - 1), x2 = stress::rand_int(x1, n - 1);
            int y1 = stress::rand_int(0, m - 1), y2 = stress::rand_int(y1, m - 1);
            check(grid, mn, mx, x1, y1, x2, y2);
        }
    }

    /// bool grids answer rectangle AND as min and OR as max, vector<bool> hands out proxies, not references
    for (int it = 0; it < 200; it++){
        int n = stress::rand_int(1, 40), m = stress::rand_int(1, 40);
        vector<vector<bool>> grid(n, vector<bool>(m));
        for (auto&& row : grid) for (auto&& x : row) x = stress::rand_int(0, 7) != 0;
        SparseTable2D<bool> mn(grid);
        SparseTable2D<bool, greater<bool>> mx(grid);
        for (int q = 0; q < 100; q++){
            int x1 = stress::rand_int(0, n - 1), x2 = stress::rand_int(x1, n - 1);
            int y1 = stress::rand_int(0, m - 1), y2 = stress::rand_int(y1, m - 1);
            check(grid, mn, mx, x1, y1, x2, y2);
        }
    }

    SparseTable2D<int> empty(vector<vector<int>>{});
    assert(empty.table.empty());

    /// The header's 500 x 500 memory figure
    auto grid = random_grid<int>(500, 500, INT_MIN, INT_MAX);
    SparseTable2D<int> mn(grid);
    size_t total = 0;
    for (auto& level : mn.table) total += level.size();
    assert(total == 3998ULL * 3998);
    for (int q = 0; q < 200; q++){
        int x1 = stress::rand_int(0, 499), x2 = stress::rand_int(x1, 499);
        int y1 = stress::rand_int(0, 499), y2 = stress::rand_int(y1, 499);
        int expected = INT_MAX;
        for (int i = x1; i <= x2; i++) for (int j = y1; j <= y2; j++) expected = min(expected, grid[i][j]);
        assert(mn.query(x1, y1, x2, y2) == expected);
    }

    return 0;
}
