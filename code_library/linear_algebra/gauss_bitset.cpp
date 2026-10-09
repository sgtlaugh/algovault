/***
 *
 * Gauss-Jordan Elimination over GF(2)
 * Linear systems, rank and inverse of matrices over GF(2), with bitset rows
 *
 * Complexity: O(min(n, m) * n * MAX / 32) for gauss and matrix_rank, O(n^2 * MAX / 32) for matrix_inverse
 * C++ bitset needs a constant size, so every row operation costs MAX / 32 whatever m is: set MAX to the largest width needed
 *
 * gauss(m, equations, res): equations in m variables, bit j of equations[i] is the coefficient of x_j, bit m the right-hand side
 *     needs m < MAX, returns -1 when there is no solution
 *     otherwise returns the number of free variables, 0 for a unique solution, and res holds one solution
 * matrix_rank(m, a): rank of the n x m matrix in bits 0 to m - 1 of the rows of a, needs m <= MAX
 * matrix_inverse(n, a, inv): for the n x n matrix in a, needs n <= MAX
 *     returns false when it is singular, otherwise inv = a^-1
 *
 * For instance, the system of linear equations modulo 2
 *
 *     (2x + y - z)   % 2 = 0  ----->  equations[0] = bits {0, 1, 1, 0}
 *     (-3x - y + 2z) % 2 = 1  ----->  equations[1] = bits {1, 1, 0, 1}
 *     (-2x + y + 2z) % 2 = 1  ----->  equations[2] = bits {0, 1, 0, 1}
 *
 * has the unique solution x = 0, y = 1, z = 1
 *
 * To make it run even faster in practice, use custom bitset class with unsigned long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const int MAX = 1024;

/// Eliminates the first cols columns of a, applying every row operation to side as well when it is given
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
vector<int> eliminate_gf2(vector<bitset<MAX>>& a, int cols, vector<bitset<MAX>>* side = nullptr){
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        while (p < n && !a[p][j]) p++;
        if (p == n) continue;

        swap(a[p], a[i]);
        if (side) swap((*side)[p], (*side)[i]);
        for (int k = 0; k < n; k++){
            if (k == i || !a[k][j]) continue;
            a[k] ^= a[i];
            if (side) (*side)[k] ^= (*side)[i];
        }

        pos[j] = i++;
    }

    return pos;
}

int gauss(int m, vector<bitset<MAX>> equations, bitset<MAX>& res){
    assert(m < MAX);
    int n = equations.size(), f_var = 0;
    vector<int> pos = eliminate_gf2(equations, m);

    res.reset();
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = 1
    for (int k = rank; k < n; k++){
        if (equations[k][m]) return -1;
    }

    for (int j = 0; j < m; j++){
        if (pos[j] == -1) f_var++;
        else res[j] = equations[pos[j]][m];
    }

    return f_var;
}

bool matrix_inverse(int n, vector<bitset<MAX>> a, vector<bitset<MAX>>& inv){
    assert((int)a.size() == n && n <= MAX);
    vector<bitset<MAX>> side(n);
    for (int i = 0; i < n; i++) side[i][i] = 1;

    vector<int> pos = eliminate_gf2(a, n, &side);
    if (n && pos[n - 1] != n - 1) return false;

    inv = side;
    return true;
}

int matrix_rank(int m, vector<bitset<MAX>> a){
    vector<int> pos = eliminate_gf2(a, m);
    return m - count(pos.begin(), pos.end(), -1);
}

vector<bitset<MAX>> rows_of(const vector<string>& bits){
    vector<bitset<MAX>> res;
    for (const string& s : bits){
        bitset<MAX> row;
        for (int j = 0; j < (int)s.size(); j++) row[j] = s[j] == '1';
        res.push_back(row);
    }
    return res;
}

int main(){
    bitset<MAX> res;
    int f_var = gauss(3, rows_of({"0110", "1101", "0101"}), res);
    assert(f_var == 0);

    assert(res[0] == 0);
    assert(res[1] == 1);
    assert(res[2] == 1);

    assert(gauss(2, rows_of({"111", "110"}), res) == -1);
    assert(gauss(3, rows_of({"1101"}), res) == 2 && res[0] == 1 && !res[1] && !res[2]);

    assert(matrix_rank(3, rows_of({"110", "011", "101"})) == 2);
    assert(matrix_rank(3, rows_of({"110", "011", "001"})) == 3);
    assert(matrix_rank(0, {}) == 0);

    /// {{1, 1}, {0, 1}} is its own inverse mod 2, and {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}} has inverse {{1, 1, 1}, {0, 1, 1}, {0, 0, 1}}
    vector<bitset<MAX>> inv;
    assert(matrix_inverse(2, rows_of({"11", "01"}), inv) && inv == rows_of({"11", "01"}));
    assert(matrix_inverse(3, rows_of({"110", "011", "001"}), inv) && inv == rows_of({"111", "011", "001"}));
    assert(!matrix_inverse(3, rows_of({"110", "011", "101"}), inv));
    assert(matrix_inverse(0, {}, inv) && inv.empty());

    return 0;
}
