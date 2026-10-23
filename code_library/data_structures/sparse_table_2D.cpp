/***
 *
 * 2D Sparse Table
 * O(1) rectangle minimum or maximum on a static grid
 *
 * Complexity: O(n m log n log m) build time and memory, O(1) per query
 *
 * SparseTable2D<T> table(grid): rectangle minimum of an n x m grid
 * SparseTable2D<T, greater<T>> table(grid): rectangle maximum
 * table.query(x1, y1, x2, y2): the answer over rows x1..x2 and columns y1..y2, 0 <= x1 <= x2 < n, 0 <= y1 <= y2 < m
 *
 * A custom Compare needs a const operator(), since query() is const
 * The grid must be rectangular, an empty grid builds but has no valid query
 * A 500 x 500 grid stores 3998^2 values, 64 MB of int
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
        assert(0 <= x1 && x1 <= x2 && x2 < n && 0 <= y1 && y1 <= y2 && y2 < m);
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
    /// Minimum and maximum over rectangles of one grid
    vector<vector<int>> grid = {
        {3, 8, 1, 9},
        {7, 2, 6, 4},
        {5, 0, 9, 3},
    };
    SparseTable2D<int> mn(grid);
    assert(mn.query(0, 0, 2, 3) == 0);
    assert(mn.query(0, 0, 1, 1) == 2);  /// rows 0..1, columns 0..1: {3, 8, 7, 2}

    SparseTable2D<int, greater<int>> mx(grid);
    assert(mx.query(0, 0, 2, 3) == 9);
    assert(mx.query(1, 0, 1, 3) == 7);

    return 0;
}
