/***
 *
 * Gauss-Jordan Elimination over GF(2)
 * Linear systems, rank and inverse of matrices over GF(2), with bitset rows
 *
 * Complexity: O(min(n, m) n W / 64) for gauss and matrix_rank, O(n^2 W / 64) for matrix_inverse
 *
 * gauss(m, equations): {solvable, x, free_variables} for equations in m < W variables, bit j of a row is the coefficient of x_j and bit m the right-hand side
 * matrix_rank(m, a): rank of the n x m matrix in bits 0 to m - 1 of the rows of a, m <= W
 * matrix_inverse(n, a): a^-1 for the n x n matrix in a, n <= W, nullopt when it is singular
 * eliminate_gf2(a, cols, side = nullptr): pos with pos[j] the pivot row of column j or -1, after reducing the first cols <= W columns of a in place
 * rows_of<W>(strings): rows from strings of '0' and '1', character j is bit j
 *
 * x is one solution, all zero when solvable is false, and free_variables is m minus the rank either way
 * eliminate_gf2 leaves reduced row echelon form with pivot rows 0, 1, ..., rank - 1 in column order, and repeats every row operation on *side
 * Rows are bitset<W>, W is deduced from the arguments and defaults to GF2_WIDTH = 4096, several widths can be used in one program
 * Every row operation costs W / 64 whatever m is, so pick the smallest W that fits: bitset<128> rows for 100 columns XOR 2 words per row operation instead of 64
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const int GF2_WIDTH = 4096;

template <size_t W>
struct GaussGF2Solution{
    bool solvable;
    bitset<W> x;
    int free_variables;
};

/// Row echelon form of the first cols columns, only rows below a pivot are cleared
template <size_t W>
vector<int> echelon_gf2(vector<bitset<W>>& a, int cols, vector<bitset<W>>* side = nullptr){
    assert(cols <= (int)W);
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        while (p < n && !a[p][j]) p++;
        if (p == n) continue;

        swap(a[p], a[i]);
        if (side) swap((*side)[p], (*side)[i]);
        for (int k = i + 1; k < n; k++){
            if (!a[k][j]) continue;
            a[k] ^= a[i];
            if (side) (*side)[k] ^= (*side)[i];
        }

        pos[j] = i++;
    }

    return pos;
}

template <size_t W>
vector<int> eliminate_gf2(vector<bitset<W>>& a, int cols, vector<bitset<W>>* side = nullptr){
    vector<int> pos = echelon_gf2(a, cols, side);
    for (int j = 0; j < cols; j++){
        if (pos[j] == -1) continue;

        int i = pos[j];
        for (int k = 0; k < i; k++){
            if (!a[k][j]) continue;
            a[k] ^= a[i];
            if (side) (*side)[k] ^= (*side)[i];
        }
    }

    return pos;
}

template <size_t W = GF2_WIDTH>
GaussGF2Solution<W> gauss(int m, vector<bitset<W>> equations){
    assert(m < (int)W);
    int n = equations.size();
    vector<int> pos = echelon_gf2(equations, m);
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = 1
    for (int k = rank; k < n; k++){
        if (equations[k][m]) return {false, {}, m - rank};
    }

    /// Back substitution from the last pivot up, x holds only the variables after j, so the AND is their parity in the row
    bitset<W> x;
    for (int j = m - 1; j >= 0; j--){
        if (pos[j] == -1) continue;

        bitset<W>& row = equations[pos[j]];
        x[j] = row[m] ^ ((row & x).count() & 1);
    }

    return {true, x, m - rank};
}

/// W is not deduced from an empty braced list, matrix_rank(0, {}), so it falls back to GF2_WIDTH
template <size_t W = GF2_WIDTH>
int matrix_rank(int m, vector<bitset<W>> a){
    vector<int> pos = echelon_gf2(a, m);
    return m - count(pos.begin(), pos.end(), -1);
}

template <size_t W = GF2_WIDTH>
vector<bitset<W>> rows_of(const vector<string>& bits){
    vector<bitset<W>> res;
    for (const string& s : bits){
        bitset<W> row;
        for (int j = 0; j < (int)s.size(); j++) row[j] = s[j] == '1';
        res.push_back(row);
    }
    return res;
}

template <size_t W>
int min_presses(int m, vector<bitset<W>> equations){
    vector<int> pos = eliminate_gf2(equations, m);
    int rank = m - count(pos.begin(), pos.end(), -1);
    for (int k = rank; k < (int)equations.size(); k++){
        if (equations[k][m]) return -1;
    }

    vector<int> free_columns;
    for (int j = 0; j < m; j++){
        if (pos[j] == -1) free_columns.push_back(j);
    }

    int best = INT_MAX;
    for (int mask = 0; mask < (1 << free_columns.size()); mask++){
        bitset<W> free_values;
        for (int b = 0; b < (int)free_columns.size(); b++){
            if (mask >> b & 1) free_values[free_columns[b]] = 1;
        }

        int presses = __builtin_popcount(mask);
        for (int i = 0; i < rank; i++) presses += equations[i][m] ^ ((equations[i] & free_values).count() & 1);
        best = min(best, presses);
    }

    return best;
}

int main(){
    int t;

    scanf("%d", &t);
    while (t--){
        int n;
        scanf("%d", &n);
        vector<string> wall(n);
        for (int i = 0; i < n; i++){
            char row[16];
            scanf("%s", row);
            wall[i] = row;
        }

        vector<bitset<256>> equations(n * n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                auto& eq = equations[i * n + j];
                eq[n * n] = wall[i][j] == 'w', eq[i * n + j] = 1;
                if (i > 0) eq[(i - 1) * n + j] = 1;
                if (i + 1 < n) eq[(i + 1) * n + j] = 1;
                if (j > 0) eq[i * n + j - 1] = 1;
                if (j + 1 < n) eq[i * n + j + 1] = 1;
            }
        }

        int presses = min_presses(n * n, equations);
        if (presses == -1) puts("inf");
        else printf("%d\n", presses);
    }
    return 0;
}
