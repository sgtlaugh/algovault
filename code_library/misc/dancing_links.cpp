/***
 *
 * Dancing Links
 * Solves exact cover with Knuth's Algorithm X: select a subset of rows covering every column exactly once
 *
 * Complexity: exponential worst case (exact cover is NP-complete), O(columns + total row length) memory
 *
 * DancingLinks(ncolumns, max_nodes): columns are 1 based, numbered 1 to ncolumns
 * max_nodes = total length of all rows presizes the node arrays, more rows still fit after a reallocation
 * add_row(r, columns): r is any int id chosen by the caller, columns must be distinct, an empty row is ignored
 * The search branches on the first column with the fewest remaining rows, stopping early at a count <= 1
 * The search is deterministic: the same rows added in the same order give the same solution
 * exact_cover restores the links before returning, so it can be called again or rows can still be added
 * Recursion depth is at most the number of selected rows, at most ncolumns
 * An 8 MB stack overflows near depth 1e5 with -O2 (near 5e4 under sanitizers), deeper covers need a larger stack
 *
 * Example:
 *   DancingLinks dlx(3);
 *   dlx.add_row(1, {1, 2});
 *   dlx.add_row(2, {3});
 *   vector<int> rows;
 *   dlx.exact_cover(rows);  // true, rows = {1, 2}
 *
 * The sudoku namespace below reduces an n x n sudoku (n a perfect square) to exact cover with 4 n^2 columns
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DancingLinks{
    int n, idx;
    vector<int> row, col, L, R, U, D, selected_rows, column_count;

    DancingLinks(int ncolumns, int max_nodes = 0) : n(ncolumns), idx(ncolumns + 1), selected_rows(ncolumns + 1), column_count(ncolumns + 1){
        for (auto* v : {&row, &col, &L, &R, &U, &D}) v->resize(n + 1 + max_nodes);

        for (int i = 0; i <= n; i++) U[i] = D[i] = i, L[i] = i - 1, R[i] = i + 1;
        L[0] = n, R[n] = 0;
    }

    void add_row(int r, const vector<int>& columns){
        int first = idx, l = columns.size();
        if (!l) return;  /// an empty row covers nothing, and linking it would corrupt the previous row
        if (idx + l > (int)row.size()){
            for (auto* v : {&row, &col, &L, &R, &U, &D}) v->resize(2 * (idx + l));
        }

        for (int c : columns){
            L[idx] = idx - 1, R[idx] = idx + 1, D[idx] = c, U[idx] = U[c];
            D[U[c]] = idx, U[c] = idx, row[idx] = r, col[idx] = c;
            column_count[c]++, idx++;
        }
        R[idx - 1] = first, L[first] = idx - 1;
    }

    /// Returns true and the selected row ids if an exact cover exists, false otherwise
    bool exact_cover(vector<int>& rows){
        int len = algorithm_x(0);
        if (len == -1) return false;

        rows.assign(selected_rows.begin(), selected_rows.begin() + len);
        return true;
    }

private:
    /// Returns the number of selected rows of the first cover found, or -1 if none exists
    int algorithm_x(int depth){
        if (R[0] == 0) return depth;

        int c = R[0];
        /// Stops at a count <= 1: the choice is forced, and a 0 count column left behind still fails the branch
        for (int i = R[0]; i != 0 && column_count[c] > 1; i = R[i]){
            if (column_count[i] < column_count[c]) c = i;
        }

        remove(c);
        int len = -1;
        for (int i = D[c]; i != c && len == -1; i = D[i]){
            selected_rows[depth] = row[i];
            for (int j = R[i]; j != i; j = R[j]) remove(col[j]);
            len = algorithm_x(depth + 1);
            for (int j = L[i]; j != i; j = L[j]) restore(col[j]);
        }

        restore(c);
        return len;
    }

    void remove(int c){
        int l = L[c], r = R[c];
        L[r] = l, R[l] = r;

        for (int i = D[c]; i != c; i = D[i]){
            for (int j = R[i]; j != i; j = R[j]){
                int u = U[j], d = D[j];
                column_count[col[j]]--;
                U[d] = u, D[u] = d;
            }
        }
    }

    void restore(int c){
        for (int i = U[c]; i != c; i = U[i]){
            for (int j = L[i]; j != i; j = L[j]){
                int u = U[j], d = D[j];
                column_count[col[j]]++;
                U[d] = j, D[u] = j;
            }
        }

        int l = L[c], r = R[c];
        L[r] = c, R[l] = c;
    }
};

namespace sudoku{
    int encode(int n, int a, int b, int c){
        return a * n * n + b * n + c + 1;
    }

    void decode(int n, int v, int& a, int& b, int& c){
        v--;
        c = v % n, v /= n;
        b = v % n, a = (v / n) % n;
    }

    /***
     * Returns false if a sudoku grid is unsolvable
     * Otherwise returns true and solves the grid in place
     * Grid must be a square and its side length must be a perfect square
     *
     * The exact cover search is exponential in the worst case
     * 9 x 9 grids, well constrained 16 x 16 grids and even empty 49 x 49 grids solve fast,
     * but sparse 16 x 16 grids (around 75% blank) can take over a minute
     *
    ***/
    bool solve(vector<vector<int>>& grid){
        int i, j, k, n = grid.size(), m = sqrt(n + 0.5);

        assert(m * m == n);
        for (i = 0; i < n; i++) assert((int)grid[i].size() == n);

        int nrows = 0;
        for (auto& r : grid) for (int x : r) nrows += x ? 1 : n;

        DancingLinks dlx(4 * n * n, 4 * nrows);  /// n * n for cells, n * n for rows, n * n for columns and n * n for boxes
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                for (k = 0; k < n; k++){
                    if (grid[i][j] == 0 || grid[i][j] == (k + 1)){
                        dlx.add_row(encode(n, i, j, k), {encode(n, 0, i, j), encode(n, 1, i, k), encode(n, 2, j, k), encode(n, 3, (i / m) * m + j / m, k)});
                    }
                }
            }
        }

        vector<int> res;
        if (!dlx.exact_cover(res)) return false;

        for (int v : res){
            decode(n, v, i, j, k);
            grid[i][j] = k + 1;
        }
        return true;
    }
}

int main(){
    DancingLinks knuth(7);
    for (auto& [r, columns] : vector<pair<int, vector<int>>>{{1, {3, 5, 6}}, {2, {1, 4, 7}}, {3, {2, 3, 6}}, {4, {1, 4}}, {5, {2, 7}}, {6, {4, 5, 7}}}){
        knuth.add_row(r, columns);
    }

    DancingLinks blocked(3);
    blocked.add_row(10, {1, 2});
    blocked.add_row(20, {2, 3});
    blocked.add_row(30, {});

    vector<int> rows;
    assert(knuth.exact_cover(rows));
    sort(rows.begin(), rows.end());
    assert((rows == vector<int>{1, 4, 5}));
    assert(!blocked.exact_cover(rows));
    assert(knuth.exact_cover(rows));
    sort(rows.begin(), rows.end());
    assert((rows == vector<int>{1, 4, 5}));

    blocked.add_row(40, {3});
    assert(blocked.exact_cover(rows));
    sort(rows.begin(), rows.end());
    assert((rows == vector<int>{10, 40}));

    DancingLinks nothing(0);
    assert(nothing.exact_cover(rows) && rows.empty());

    DancingLinks uncoverable(2);
    uncoverable.add_row(7, {1});
    assert(!uncoverable.exact_cover(rows));

    vector<vector<int>> grid = {
        {0, 0, 0, 0, 6, 9, 8, 3, 0},
        {9, 8, 0, 0, 0, 0, 0, 7, 6},
        {6, 0, 0, 0, 3, 8, 0, 5, 1},
        {2, 0, 5, 0, 8, 1, 0, 9, 0},
        {0, 6, 0, 0, 0, 0, 0, 8, 0},
        {0, 9, 0, 3, 7, 0, 6, 0, 2},
        {3, 4, 0, 8, 5, 0, 0, 0, 9},
        {7, 2, 0, 0, 0, 0, 0, 6, 8},
        {0, 5, 6, 9, 2, 0, 0, 0, 0},
    };

    assert(sudoku::solve(grid));

    const vector<vector<int>> solved_grid = {
        {5, 1, 2, 7, 6, 9, 8, 3, 4},
        {9, 8, 3, 5, 1, 4, 2, 7, 6},
        {6, 7, 4, 2, 3, 8, 9, 5, 1},
        {2, 3, 5, 6, 8, 1, 4, 9, 7},
        {1, 6, 7, 4, 9, 2, 3, 8, 5},
        {4, 9, 8, 3, 7, 5, 6, 1, 2},
        {3, 4, 1, 8, 5, 6, 7, 2, 9},
        {7, 2, 9, 1, 4, 3, 5, 6, 8},
        {8, 5, 6, 9, 2, 7, 1, 4, 3},
    };
    assert(grid == solved_grid);

    grid = {
        {0, 0, 0, 0, 6, 9, 8, 3, 0},
        {9, 8, 0, 0, 0, 0, 0, 7, 6},
        {6, 0, 0, 0, 3, 8, 0, 5, 1},
        {2, 0, 5, 4, 8, 1, 0, 9, 0},
        {0, 6, 0, 0, 0, 0, 0, 8, 0},
        {0, 9, 0, 3, 7, 0, 6, 0, 2},
        {3, 4, 0, 8, 5, 0, 0, 0, 9},
        {7, 2, 0, 0, 0, 0, 0, 6, 8},
        {0, 5, 6, 9, 2, 0, 0, 0, 0},
    };
    assert(!sudoku::solve(grid));

    int n = 25;
    auto large_grid = vector<vector<int>>(n, vector<int>(n, 0));
    assert(sudoku::solve(large_grid));

    const vector<vector<int>> solved_large_grid = {
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25},
        {11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 16, 17, 18, 19, 20, 22, 23, 21, 24, 25, 6, 7, 8, 9, 10},
        {6, 7, 8, 9, 10, 24, 25, 22, 21, 23, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
        {23, 22, 21, 25, 24, 16, 17, 18, 19, 20, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 11, 12, 13, 14, 15},
        {16, 17, 18, 19, 20, 11, 12, 13, 14, 15, 24, 21, 23, 22, 25, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5},
        {3, 1, 5, 2, 14, 19, 16, 20, 7, 11, 18, 15, 21, 24, 23, 8, 10, 25, 22, 12, 4, 9, 17, 13, 6},
        {4, 8, 10, 6, 18, 3, 1, 23, 24, 2, 13, 22, 25, 12, 14, 9, 16, 19, 11, 17, 5, 15, 20, 21, 7},
        {24, 15, 25, 7, 12, 4, 5, 9, 13, 6, 3, 1, 16, 17, 2, 18, 20, 14, 21, 23, 8, 10, 22, 11, 19},
        {9, 19, 11, 20, 21, 8, 22, 25, 17, 14, 4, 5, 10, 7, 6, 3, 1, 15, 13, 2, 18, 24, 16, 12, 23},
        {13, 23, 17, 22, 16, 10, 15, 21, 18, 12, 8, 9, 20, 11, 19, 4, 5, 24, 7, 6, 3, 1, 14, 25, 2},
        {2, 5, 1, 15, 3, 14, 13, 17, 25, 8, 23, 16, 19, 10, 18, 12, 24, 11, 20, 9, 7, 4, 21, 6, 22},
        {7, 16, 6, 11, 4, 2, 18, 1, 5, 3, 12, 24, 17, 8, 22, 14, 21, 23, 10, 13, 20, 25, 19, 15, 9},
        {12, 25, 23, 21, 8, 7, 4, 6, 11, 9, 2, 20, 1, 5, 3, 17, 15, 22, 18, 19, 14, 13, 24, 10, 16},
        {14, 20, 19, 24, 9, 22, 23, 10, 16, 21, 7, 4, 6, 15, 13, 2, 25, 1, 5, 3, 12, 8, 11, 17, 18},
        {17, 10, 22, 18, 13, 12, 24, 15, 20, 19, 14, 25, 11, 21, 9, 7, 4, 6, 8, 16, 2, 23, 1, 5, 3},
        {19, 6, 4, 1, 2, 25, 21, 11, 12, 24, 20, 18, 9, 23, 8, 15, 14, 16, 17, 22, 10, 5, 7, 3, 13},
        {5, 13, 9, 3, 7, 17, 8, 2, 1, 4, 22, 19, 24, 25, 11, 10, 18, 20, 23, 21, 15, 6, 12, 16, 14},
        {15, 18, 20, 23, 17, 5, 6, 7, 3, 16, 10, 14, 2, 1, 4, 19, 13, 9, 12, 24, 22, 21, 25, 8, 11},
        {22, 14, 12, 16, 11, 18, 10, 19, 15, 13, 5, 6, 7, 3, 21, 25, 8, 2, 1, 4, 17, 20, 9, 23, 24},
        {10, 21, 24, 8, 25, 20, 9, 14, 23, 22, 15, 13, 12, 16, 17, 5, 6, 7, 3, 11, 19, 18, 2, 1, 4},
        {18, 3, 2, 5, 1, 15, 20, 4, 8, 25, 17, 10, 22, 13, 24, 23, 11, 12, 16, 14, 9, 19, 6, 7, 21},
        {25, 4, 7, 13, 6, 9, 3, 16, 2, 1, 19, 23, 14, 20, 12, 21, 22, 10, 15, 8, 24, 11, 5, 18, 17},
        {20, 11, 15, 12, 19, 23, 14, 5, 6, 7, 21, 3, 4, 2, 1, 24, 9, 17, 25, 18, 13, 16, 10, 22, 8},
        {8, 9, 16, 17, 22, 21, 19, 24, 10, 18, 25, 11, 5, 6, 7, 13, 3, 4, 2, 1, 23, 14, 15, 20, 12},
        {21, 24, 14, 10, 23, 13, 11, 12, 22, 17, 9, 8, 15, 18, 16, 20, 19, 5, 6, 7, 25, 3, 4, 2, 1},
    };
    assert(large_grid == solved_large_grid);

    clock_t start = clock();

    n = 49;
    auto very_large_grid = vector<vector<int>>(n, vector<int>(n, 0));
    assert(sudoku::solve(very_large_grid));

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));  /// 0.31600 seconds

    return 0;
}
