/***
 *
 * Dancing Links
 * Solves exact cover with Knuth's Algorithm X: select a subset of rows covering every column exactly once
 *
 * Complexity: exponential worst case (exact cover is NP-complete), O(columns + total row length) memory
 *
 * DancingLinks dlx(ncolumns, max_nodes = 0): columns 0 to ncolumns - 1, max_nodes = total row length presizes arrays
 * dlx.add_row(r, columns): adds a row with any int id r over distinct columns, an empty row is ignored
 * dlx.count(limit = LLONG_MAX): number of exact covers, or limit if there are more, searching no further, limit >= 0
 * dlx.exact_cover(): the selected row ids, nullopt if no cover exists, an empty cover for 0 columns
 *
 * Notes:
 *   - Rows beyond max_nodes still fit after a reallocation
 *   - Deterministic: the same rows added in the same order give the same cover
 *   - count(2) == 1 tells a unique cover, like a sudoku with one solution, without searching for all of them
 *   - count and exact_cover restore the links, so they can be called again or rows can still be added
 *   - Recursion depth is at most ncolumns, an 8 MB stack overflows near 5e4 (4e4 under sanitizers)
 *   - sudoku::solution(grid) below reduces an n x n sudoku (n a perfect square) to exact cover with 4 n^2 columns
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DancingLinks{
    int n, idx, cover_size = 0;
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

        for (int column : columns){
            assert(0 <= column && column < n);
            int c = column + 1;  /// node 0 is the header of the column list
            L[idx] = idx - 1, R[idx] = idx + 1, D[idx] = c, U[idx] = U[c];
            D[U[c]] = idx, U[c] = idx, row[idx] = r, col[idx] = c;
            column_count[c]++, idx++;
        }
        R[idx - 1] = first, L[first] = idx - 1;
    }

    long long count(long long limit = LLONG_MAX){
        assert(limit >= 0);
        return limit ? algorithm_x(0, limit) : 0;
    }

    optional<vector<int>> exact_cover(){
        if (!algorithm_x(0, 1)) return nullopt;
        return vector<int>(selected_rows.begin(), selected_rows.begin() + cover_size);
    }

    /// Counts covers up to limit >= 1, the last one found stays in selected_rows[0, cover_size)
    long long algorithm_x(int depth, long long limit){
        if (R[0] == 0){
            cover_size = depth;
            return 1;
        }

        int c = R[0];
        /// Stops at a count <= 1: the choice is forced, and a 0 count column left behind still fails the branch
        for (int i = R[0]; i != 0 && column_count[c] > 1; i = R[i]){
            if (column_count[i] < column_count[c]) c = i;
        }

        remove(c);
        long long found = 0;
        for (int i = D[c]; i != c && found < limit; i = D[i]){
            selected_rows[depth] = row[i];
            for (int j = R[i]; j != i; j = R[j]) remove(col[j]);
            found += algorithm_x(depth + 1, limit - found);
            for (int j = L[i]; j != i; j = L[j]) restore(col[j]);
        }

        restore(c);
        return found;
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
        return a * n * n + b * n + c;
    }

    void decode(int n, int v, int& a, int& b, int& c){
        c = v % n, v /= n;
        b = v % n, a = (v / n) % n;
    }

    /***
     * Returns the solved grid, nullopt if the puzzle is unsolvable; the puzzle itself is not changed
     * Zeros are blanks, the grid must be n x n with n a perfect square
     *
     * The exact cover search is exponential in the worst case
     * 9 x 9 grids, well constrained 16 x 16 grids and even empty 49 x 49 grids solve fast,
     * but sparse 16 x 16 grids (around 75% blank) can take over a minute
     *
    ***/
    optional<vector<vector<int>>> solution(const vector<vector<int>>& grid){
        int i, j, k, n = grid.size(), m = sqrt(n + 0.5);

        assert(m * m == n);
        for (i = 0; i < n; i++) assert((int)grid[i].size() == n);

        int nrows = 0;
        for (auto& r : grid) for (int x : r) nrows += x ? 1 : n;

        DancingLinks dlx(4 * n * n, 4 * nrows);  /// n * n columns each for cells, rows, columns and boxes
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                for (k = 0; k < n; k++){
                    if (grid[i][j] == 0 || grid[i][j] == (k + 1)){
                        dlx.add_row(encode(n, i, j, k), {encode(n, 0, i, j), encode(n, 1, i, k), encode(n, 2, j, k), encode(n, 3, (i / m) * m + j / m, k)});
                    }
                }
            }
        }

        auto res = dlx.exact_cover();
        if (!res) return nullopt;

        vector<vector<int>> solved = grid;
        for (int v : *res){
            decode(n, v, i, j, k);
            solved[i][j] = k + 1;
        }

        return solved;
    }
}

int main(){
    /// Knuth's example: rows 1, 4 and 5 cover the 7 columns, each exactly once
    DancingLinks knuth(7);
    knuth.add_row(1, {2, 4, 5});
    knuth.add_row(2, {0, 3, 6});
    knuth.add_row(3, {1, 2, 5});
    knuth.add_row(4, {0, 3});
    knuth.add_row(5, {1, 6});
    knuth.add_row(6, {3, 4, 6});

    auto rows = knuth.exact_cover();
    assert(rows);
    sort(rows->begin(), rows->end());
    assert((*rows == vector<int>{1, 4, 5}));

    /// Column 0 forces row 10, then column 2 has no row until 40 is added
    DancingLinks blocked(3);
    blocked.add_row(10, {0, 1});
    blocked.add_row(20, {1, 2});
    assert(!blocked.exact_cover());

    blocked.add_row(40, {2});
    rows = blocked.exact_cover();
    assert(rows);
    assert((*rows == vector<int>{10, 40}));

    /// Two covers: {1} alone, or {2, 3}
    DancingLinks twice(2);
    twice.add_row(1, {0, 1});
    twice.add_row(2, {0});
    twice.add_row(3, {1});
    assert(twice.count() == 2);
    assert(twice.count(1) == 1);

    const vector<vector<int>> grid = {
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

    auto solved = sudoku::solution(grid);
    assert(solved);

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
    assert((*solved == solved_grid));

    return 0;
}
