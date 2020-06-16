/***

 * Dancing Links data structure to solve exact cover problems using Algorithm X

 * There are some constraints as columns and a number of rows (both 1 based)
 * Each row satisfies some constraints
 * Objective is to select a subset of rows so that each constraint is satisfied exactly once
 * Don't forget to initialize first by calling init()

***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

/// Define MAX limits appropriately, set to large values for safety
#define MAXR 2000010
#define MAXC 2000010
#define MAXNODE 2000010

using namespace std;

namespace dlx{
    int row[MAXNODE], col[MAXNODE];
    int L[MAXNODE], R[MAXNODE], U[MAXNODE], D[MAXNODE];
    int n, idx, len, selected_rows[MAXR], column_count[MAXC];

    void init(int ncolumn){ /// initialize first with total number of columns (1 based)
        memset(column_count, 0, sizeof(column_count));

        n = ncolumn, idx = n + 1;
        for (int i = 0; i <= n; i++) U[i] = D[i] = i, L[i] = i - 1, R[i] = i + 1;
        L[0] = n, R[n] = 0;
    }

    /// r = index of row (1 based)
    /// the vector columns contain the columns which are satisfied with this row
    inline void addrow(int r, vector <int>& columns){
        int i, c, l = columns.size(), first = idx;

        for (i = 0; i < l; i++){
            c = columns[i];
            L[idx] = idx - 1, R[idx] = idx + 1, D[idx] = c, U[idx] = U[c];
            D[U[c]] = idx, U[c] = idx, row[idx] = r, col[idx] = c;
            column_count[c]++, idx++; /// column_count[c] is the number of rows which satisfies constraint column c
        }
        R[idx - 1] = first, L[first] = idx - 1;
    }

    /// Removes column c from the structure
    inline void remove(int c){
        L[R[c]] = L[c], R[L[c]] = R[c];

        for (int i = D[c]; i != c; i = D[i]){
            for (int j = R[i]; j != i; j = R[j]){
                column_count[col[j]]--;
                U[D[j]] = U[j], D[U[j]] = D[j];
            }
        }
    }

    /// Restores the position of column c in the structure
    inline void restore(int c){
        for (int i = U[c]; i != c; i = U[i]){
            for (int j = L[i]; j != i; j = L[j]){
                column_count[col[j]]++;
                U[D[j]] = j, D[U[j]] = j;
            }
        }

        L[R[c]] = c, R[L[c]] = c;
    }

    /// Recursively enumerate to solve exact cover
    bool algorithmX(int depth){
        if(R[0] == 0){
            len = depth;
            return true;
        }

        int i, j, c = R[0];
        for (i = R[0]; i != 0; i = R[i]){ /// Select a column deterministically
            if(column_count[i] < column_count[c]) c = i;
        }

        remove(c);
        bool flag = false;
        for (i = D[c]; i != c && !flag; i = D[i]){
            selected_rows[depth] = row[i];
            for (j = R[i]; j != i; j = R[j]) remove(col[j]);
            flag |= algorithmX(depth + 1); /// Perhaps select rows non-deterministically here with random_shuffle for optimizations?
            for (j = L[i]; j != i; j = L[j]) restore(col[j]);
        }

        restore(c);
        return flag;
    }

    /// Returns the subset of rows satisfying exact cover, false otherwise
    bool exact_cover(vector<int>& rows){
        if(!algorithmX(0)) return false;
        rows = vector<int>(selected_rows, selected_rows + len);
        return true;
    }
}

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
     * Returns false if a sudoku grid is solvable
     * Otherwise returns true and solves the grid in place
     * Grid must be a square and its side length must be a perfect square
     *
     * Algorithm should be fast enough for any 16 x 16 grid or lower
     * Although certain configurations of higher grid sizes are definitely solvable
     *
    ***/

    bool solve(vector <vector<int>>& grid){
        int i, j, k, l, n = grid.size(), m = sqrt(n + 0.5);

        assert(m * m == n);
        for (i = 0; i < n; i++) assert((int)grid[i].size() == n);

        dlx::init(4 * n * n); /// n * n for cells, n * n for rows, n * n for columns and n * n for boxes
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                for (k = 0; k < n; k++){
                    if (grid[i][j] == 0 || grid[i][j] == (k + 1)){
                        vector <int> columns;
                        columns.push_back(encode(n, 0, i, j));
                        columns.push_back(encode(n, 1, i, k));
                        columns.push_back(encode(n, 2, j, k));
                        columns.push_back(encode(n, 3, (i / m) * m + j / m, k));
                        dlx::addrow(encode(n, i, j, k), columns);
                    }
                }
            }
        }

        vector<int> res;
        if (!dlx::exact_cover(res)) return false;

        for (l = 0; l < (int)res.size(); l++){
            decode(n, res[l], i, j, k);
            grid[i][j] = k + 1;
        }
        return true;
    }
}

int main(){
    vector <vector<int>> grid = {
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
