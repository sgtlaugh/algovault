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

    T best(const T& x, const T& y) const{
        return compare(y, x) ? y : x;
    }

    /// By value, vector<bool>::operator[] returns a temporary that a reference would outlive
    T block(int a, int b, int i, int j) const{
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

    assert(mn.query(0, 0, 2, 3) == 0);      /// the whole grid
    assert(mx.query(0, 0, 2, 3) == 9);
    assert(mn.query(0, 0, 1, 1) == 2);      /// rows 0..1, columns 0..1: {3, 8, 7, 2}
    assert(mx.query(1, 0, 1, 3) == 7);      /// row 1 alone
    assert(mx.query(2, 1, 2, 1) == 0);      /// a single cell
    return 0;
}
