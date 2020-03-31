#include <bits/stdtr1c++.h>

#define MAXR 100010
#define MAXC 100010
#define MAXNODE 100010
/// Define MAX limits appropriately, set to large values for safety

using namespace std;

/***

Dancing Links data structure to solve exact cover problems using Algorithm X

There are some constraints as columns and a number of rows
Each row satisfies some constraints
Objective is to select a subset of rows so that each constraint is satisfied exactly once
Don't forget to initialize first by calling init()

***/

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
            if(column_count[i] < column_count[c]) c = i; /// if multiple columns exist, choose the one with minimum count
        }

        remove(c);
        bool flag = false;
        for (i = D[c]; i != c && !flag; i = D[i]){ /// While solution is not found
            selected_rows[depth] = row[i];
            for (j = R[i]; j != i; j = R[j]) remove(col[j]);
            flag |= algorithmX(depth + 1); /// May be select rows non-deterministically here with random_shuffle?
            for (j = L[i]; j != i; j = L[j]) restore(col[j]);
        }

        restore(c);
        return flag;
    }

    bool exact_cover(vector<int>& rows){ /// Returns the subset of rows satisfying exact cover, false otherwise
        rows.clear();
        if(!algorithmX(0)) return false;
        for(int i = 0; i < len; i++) rows.push_back(selected_rows[i]);
        return true;
    }
}

namespace sudoku{
    inline int encode(int n, int a, int b, int c){ /// Encode information
        return (a * n * n) + (b * n) + c + 1;
    }

    inline void decode(int n, int v, int& a, int& b, int& c){ /// Decode information
        v--;
        c = v % n, v /= n;
        b = v % n, a = (v / n) % n;
    }

    bool solve(int n, int ar[25][25]){
        int i, j, k, l;
        int m = sqrt(n + 0.5); /// m * m sub-grid
        int ncolumn = n * n * 4; /// n * n for cells, n * n for rows, n * n for columns and n * n for boxes

        dlx::init(ncolumn);
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                for (k = 0; k < n; k++){
                    if (ar[i][j] == 0 || ar[i][j] == (k + 1)){
                        vector <int> columns;
                        columns.push_back(encode(n, 0, i, j)); /// Each cell can contain only 1 value
                        columns.push_back(encode(n, 1, i, k)); /// Row[i] can contain only ar[i][j], if ar[i][j] != 0, otherwise it can contain any value
                        columns.push_back(encode(n, 2, j, k)); /// Column[j] can contain only ar[i][j], if ar[i][j] != 0, otherwise it can contain any value
                        columns.push_back(encode(n, 3, (i / m) * m + j / m, k)); /// Similarly for box numbered i / m * m + j / m

                        /// Current row selection, place number k in the cell[i][j] and doing so will cover all columns in the columns vector
                        int cur_row = encode(n, i, j, k);
                        dlx::addrow(cur_row, columns);
                    }
                }
            }
        }

        vector<int> res;
        if (!dlx::exact_cover(res)) return false;

        for (l = 0; l < (int)res.size(); l++){
            decode(n, res[l], i, j, k);
            ar[i][j] = k + 1;
        }
        return true;
    }
}

int main(){

}
