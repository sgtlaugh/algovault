/***
 *
 * Dancing Links
 * Solves exact cover with Knuth's Algorithm X: select a subset of rows covering every column exactly once
 *
 * Complexity: exponential worst case (exact cover is NP-complete), O(columns + total row length) memory
 *
 * DancingLinks(ncolumns, max_nodes): columns are 0 based, numbered 0 to ncolumns - 1
 * max_nodes = total length of all rows presizes the node arrays, more rows still fit after a reallocation
 * add_row(r, columns): r is any int id chosen by the caller, columns must be distinct, an empty row is ignored
 * The search branches on the first column with the fewest remaining rows, stopping early at a count <= 1
 * The search is deterministic: the same rows added in the same order give the same solution
 * exact_cover() returns the selected row ids, nullopt if no cover exists (0 columns give an empty cover)
 * exact_cover restores the links before returning, so it can be called again or rows can still be added
 * Recursion depth is at most the number of selected rows, at most ncolumns
 * An 8 MB stack overflows near depth 5e4 with -O2 and as early as 4e4 under sanitizers, deeper covers need a larger stack
 *
 * Example:
 *   DancingLinks dlx(3);
 *   dlx.add_row(1, {0, 1});
 *   dlx.add_row(2, {2});
 *   auto rows = dlx.exact_cover();  // {1, 2}
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

        for (int column : columns){
            assert(0 <= column && column < n);
            int c = column + 1;  /// node 0 is the header of the column list
            L[idx] = idx - 1, R[idx] = idx + 1, D[idx] = c, U[idx] = U[c];
            D[U[c]] = idx, U[c] = idx, row[idx] = r, col[idx] = c;
            column_count[c]++, idx++;
        }
        R[idx - 1] = first, L[first] = idx - 1;
    }

    optional<vector<int>> exact_cover(){
        int len = algorithm_x(0);
        if (len == -1) return nullopt;
        return vector<int>(selected_rows.begin(), selected_rows.begin() + len);
    }

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
        return a * n * n + b * n + c;
    }

    void decode(int n, int v, int& a, int& b, int& c){
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

        auto res = dlx.exact_cover();
        if (!res) return false;

        for (int v : *res){
            decode(n, v, i, j, k);
            grid[i][j] = k + 1;
        }
        return true;
    }
}

/// Counts every exact cover with the template's links: algorithm_x stops at the first cover, the problem wants all of them
int count_covers(DancingLinks& dlx){
    if (dlx.R[0] == 0) return 1;

    int c = dlx.R[0];
    for (int i = dlx.R[0]; i != 0; i = dlx.R[i]){
        if (dlx.column_count[i] < dlx.column_count[c]) c = i;
    }

    int res = 0;
    dlx.remove(c);
    for (int i = dlx.D[c]; i != c; i = dlx.D[i]){
        for (int j = dlx.R[i]; j != i; j = dlx.R[j]) dlx.remove(dlx.col[j]);
        res += count_covers(dlx);
        for (int j = dlx.L[i]; j != i; j = dlx.L[j]) dlx.restore(dlx.col[j]);
    }
    dlx.restore(c);
    return res;
}

int main(){
    int t;
    scanf("%d", &t);
    while (t--){
        char grid[9][16];
        for (int i = 0; i < 9; i++) scanf("%s", grid[i]);

        auto place = [](vector<int>& columns, int i, int j, int k){
            columns.push_back(sudoku::encode(9, 0, i, j));
            columns.push_back(sudoku::encode(9, 1, i, k));
            columns.push_back(sudoku::encode(9, 2, j, k));
            columns.push_back(sudoku::encode(9, 3, i / 3 * 3 + j / 3, k));
        };

        DancingLinks dlx(4 * 81, 4 * 81 * 9);
        vector<vector<pair<int, int>>> cells_of(26);
        for (int i = 0; i < 9; i++){
            for (int j = 0; j < 9; j++){
                char ch = grid[i][j];
                if (islower(ch) && ch != 'e' && ch != 'o'){
                    cells_of[ch - 'a'].push_back({i, j});
                    continue;
                }
                for (int k = 0; k < 9; k++){
                    int digit = k + 1;
                    if ((ch == 'e' && digit % 2) || (ch == 'o' && digit % 2 == 0) || (isdigit(ch) && ch != '0' && ch - '0' != digit)) continue;
                    vector<int> columns;
                    place(columns, i, j, k);
                    dlx.add_row(0, columns);
                }
            }
        }

        /// One row fills every cell of a variable with the same digit, so the cover keeps them equal
        for (auto& cells : cells_of){
            if (cells.empty()) continue;
            for (int k = 0; k < 9; k++){
                vector<int> columns;
                for (auto [i, j] : cells) place(columns, i, j, k);
                vector<int> sorted = columns;
                sort(sorted.begin(), sorted.end());
                /// Two cells of the variable share a row, column or box, so the same digit twice is impossible
                if (adjacent_find(sorted.begin(), sorted.end()) == sorted.end()) dlx.add_row(0, columns);
            }
        }
        printf("%d\n", count_covers(dlx));
    }
    return 0;
}
