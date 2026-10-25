/***
 *
 * Gauss-Jordan Elimination over GF(2)
 * Linear systems, rank and inverse of matrices over GF(2), with bitset rows
 *
 * Complexity: O(min(n, m) n W / 64) for gauss and matrix_rank, O(n^2 W / 64) for matrix_inverse
 *
 * eliminate_gf2(a, cols, side = nullptr): pos[j] the pivot row of column j or -1, reducing the first cols <= W columns
 * gauss(m, equations): {solvable, x, free_variables} for equations in m < W variables, bit m of a row is its right side
 * matrix_inverse(n, a): a^-1 for the n x n matrix in a, n <= W, nullopt when it is singular
 * matrix_rank(m, a): rank of the n x m matrix in bits 0 to m - 1 of the rows of a, m <= W
 * rows_of<W>(strings): rows from strings of '0' and '1', character j is bit j
 *
 * Notes:
 *   - x is one solution, all zero when solvable is false, and free_variables is m minus the rank either way
 *   - eliminate_gf2 leaves reduced row echelon form with pivot rows 0, 1, ..., rank - 1 in column order
 *   - eliminate_gf2 repeats every row operation on *side
 *   - Rows are bitset<W>, W is deduced from the arguments, default GF2_WIDTH = 4096, widths may mix in one program
 *   - A row operation costs W / 64 words whatever m is: pick the smallest W that fits, like bitset<128> for 100 columns
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

    /// Back substitution from the last pivot up: x holds only variables after j, so the AND is their parity in the row
    bitset<W> x;
    for (int j = m - 1; j >= 0; j--){
        if (pos[j] == -1) continue;

        bitset<W>& row = equations[pos[j]];
        x[j] = row[m] ^ ((row & x).count() & 1);
    }

    return {true, x, m - rank};
}

template <size_t W = GF2_WIDTH>
optional<vector<bitset<W>>> matrix_inverse(int n, vector<bitset<W>> a){
    /// Checked before the identity is written, side[i][i] with i >= W is past the end of the row
    assert((int)a.size() == n && n <= (int)W);
    vector<bitset<W>> side(n);
    for (int i = 0; i < n; i++) side[i][i] = 1;

    vector<int> pos = eliminate_gf2(a, n, &side);
    if (n && pos[n - 1] != n - 1) return nullopt;
    return side;
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

int main(){
    /// 2x + y - z = 0, -3x - y + 2z = 1, -2x + y + 2z = 1 are y + z = 0, x + y = 1, y = 1 mod 2: x = 0, y = 1, z = 1
    auto [solvable, x, free_variables] = gauss(3, rows_of({"0110", "1101", "0101"}));
    assert(solvable);
    assert(free_variables == 0);
    assert(!x[0]);
    assert(x[1]);
    assert(x[2]);

    /// x + y = 1 leaves y and z free at 0, so x = 1
    GaussGF2Solution<GF2_WIDTH> line = gauss(3, rows_of({"1101"}));
    assert(line.solvable);
    assert(line.free_variables == 2);
    assert(line.x[0]);
    assert(!line.x[1]);
    assert(!line.x[2]);

    /// x + y = 1 and x + y = 0 contradict each other
    assert(!gauss(2, rows_of({"111", "110"})).solvable);

    /// The third row is the sum of the first two
    assert(matrix_rank(3, rows_of({"110", "011", "101"})) == 2);

    /// {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}} has the inverse {{1, 1, 1}, {0, 1, 1}, {0, 0, 1}}
    optional<vector<bitset<GF2_WIDTH>>> inv = matrix_inverse(3, rows_of({"110", "011", "001"}));
    assert(inv.has_value());
    assert(*inv == rows_of({"111", "011", "001"}));

    assert(!matrix_inverse(3, rows_of({"110", "011", "101"})).has_value());

    /// x1 + x2 and x0 + x1: the pivot of column 0 moves up, then row 0 clears column 1 into x0 + x2
    vector<bitset<GF2_WIDTH>> rows = rows_of({"011", "110"});
    vector<int> pos = eliminate_gf2(rows, 3);
    assert((pos == vector<int>{0, 1, -1}));
    assert(rows == rows_of({"101", "011"}));

    /// The smallest width that fits: 3 variables plus the right-hand side in bitset<4>
    GaussGF2Solution<4> narrow = gauss(3, rows_of<4>({"0110", "1101", "0101"}));
    assert(narrow.x == bitset<4>("0110"));

    return 0;
}
