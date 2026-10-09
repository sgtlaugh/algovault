/***
 *
 * 2D Sparse Table
 * Rectangle minimum or maximum on a static n x m grid
 *
 * Complexity: O(n m log n log m) build time and memory, O(1) per query
 *
 * table[a][b] holds the answer for every 2^a x 2^b block, a query covers its rectangle with four overlapping blocks
 * Overlaps are only correct for idempotent operations, so this answers min and max, never sum or xor
 *
 * SparseTable2D<T> answers minimum, SparseTable2D<T, greater<T>> answers maximum
 * A custom Compare must have a const operator(), since query() is const
 * query(x1, y1, x2, y2) is the inclusive rectangle rows [x1, x2], columns [y1, y2], 0-indexed
 * The grid must be rectangular, an empty grid builds but has no valid query
 *
 * Stores (sum of n - 2^a + 1 over a) * (sum of m - 2^b + 1 over b) values
 * A 500 x 500 grid stores 3998^2 values, 64 MB of int
 *
 * Example:
 *   SparseTable2D<int, greater<int>> table(grid);
 *   int best = table.query(1, 0, 3, 2);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Compare = less<T>>
struct SparseTable2D{
    int n, m, lgm;
    Compare compare;
    /// table[a * lgm + b] is row-major with m - 2^b + 1 columns
    vector<vector<T>> table;

    SparseTable2D(const vector<vector<T>>& grid) : n(grid.size()), m(grid.empty() ? 0 : grid[0].size()){
        int lgn = n ? __lg(n) + 1 : 0;
        lgm = m ? __lg(m) + 1 : 0;
        table.resize(lgn * lgm);

        for (int a = 0; a < lgn; a++){
            for (int b = 0; b < lgm; b++){
                int rows = n - (1 << a) + 1, cols = m - (1 << b) + 1;
                vector<T>& cur = table[a * lgm + b];
                cur.reserve((size_t)rows * cols);

                for (int i = 0; i < rows; i++){
                    for (int j = 0; j < cols; j++){
                        if (a == 0 && b == 0) cur.push_back(grid[i][j]);
                        else if (a == 0) cur.push_back(best(block(0, b - 1, i, j), block(0, b - 1, i, j + (1 << (b - 1)))));
                        else cur.push_back(best(block(a - 1, b, i, j), block(a - 1, b, i + (1 << (a - 1)), j)));
                    }
                }
            }
        }
    }

    T query(int x1, int y1, int x2, int y2) const{
        int a = __lg(x2 - x1 + 1), b = __lg(y2 - y1 + 1);
        int x3 = x2 - (1 << a) + 1, y3 = y2 - (1 << b) + 1;
        return best(best(block(a, b, x1, y1), block(a, b, x1, y3)), best(block(a, b, x3, y1), block(a, b, x3, y3)));
    }

    const T& best(const T& x, const T& y) const{
        return compare(y, x) ? y : x;
    }

    const T& block(int a, int b, int i, int j) const{
        return table[a * lgm + b][(size_t)i * (m - (1 << b) + 1) + j];
    }
};

int main(){
    vector<vector<int>> grid = {
        {3, 8, 1, 9},
        {7, 2, 6, 4},
        {5, 0, 9, 3},
    };
    SparseTable2D<int> mn(grid);
    SparseTable2D<int, greater<int>> mx(grid);

    assert(mn.query(0, 0, 2, 3) == 0 && mx.query(0, 0, 2, 3) == 9);
    assert(mn.query(0, 0, 1, 1) == 2 && mx.query(0, 0, 1, 1) == 8);
    assert(mn.query(0, 2, 0, 3) == 1 && mx.query(0, 2, 0, 3) == 9);
    assert(mn.query(1, 2, 2, 3) == 3 && mx.query(1, 2, 2, 3) == 9);
    assert(mn.query(1, 0, 1, 3) == 2 && mx.query(1, 0, 1, 3) == 7);
    assert(mn.query(0, 0, 2, 0) == 3 && mx.query(0, 0, 2, 0) == 7);
    assert(mn.query(1, 1, 2, 3) == 0 && mx.query(0, 1, 1, 2) == 8);
    assert(mn.query(2, 1, 2, 1) == 0 && mx.query(2, 1, 2, 1) == 0);

    SparseTable2D<long long> single(vector<vector<long long>>{{LLONG_MIN}});
    assert(single.query(0, 0, 0, 0) == LLONG_MIN);

    SparseTable2D<int> empty(vector<vector<int>>{});
    assert(empty.table.empty());

    vector<vector<int>> column = {{4}, {INT_MAX}, {-2}, {INT_MIN}, {7}};
    SparseTable2D<int> col_min(column);
    SparseTable2D<int, greater<int>> col_max(column);
    assert(col_min.query(0, 0, 4, 0) == INT_MIN && col_max.query(0, 0, 4, 0) == INT_MAX);
    assert(col_min.query(0, 0, 2, 0) == -2 && col_max.query(2, 0, 4, 0) == 7);

    mt19937 rng(20020523);
    int n = 37, m = 50;
    vector<vector<int>> big(n, vector<int>(m));
    for (auto& row : big) for (auto& x : row) x = rng() % 1000;
    SparseTable2D<int, greater<int>> big_max(big);

    for (int q = 0; q < 2000; q++){
        int x1 = rng() % n, x2 = rng() % n, y1 = rng() % m, y2 = rng() % m;
        if (x1 > x2) swap(x1, x2);
        if (y1 > y2) swap(y1, y2);

        int expected = INT_MIN;
        for (int i = x1; i <= x2; i++) for (int j = y1; j <= y2; j++) expected = max(expected, big[i][j]);
        assert(big_max.query(x1, y1, x2, y2) == expected);
    }

    return 0;
}
