/***
 *
 * Gauss-Jordan Elimination modulo a prime
 * Linear systems, solution spaces, rank and inverse of matrices modulo a prime
 *
 * Complexity: O(min(n, m) n m + m^2) for gauss, O(min(n, m) n m) for matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, mod): {solvable, x, basis} for n >= 1 equations in m variables, the right side last
 * matrix_inverse(a, mod): a^-1 for a square a, nullopt when a is singular
 * matrix_rank(a, mod): rank of any n x m matrix
 *
 * Notes:
 *   - x is one solution, empty when solvable is false
 *   - basis spans {x : Ax = 0}, every solution is x plus exactly one combination of the basis vectors
 *   - basis is filled even when there is no solution, it depends only on the coefficients
 *   - mod must be a prime below 2^31, entries may be any int and results are in [0, mod)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GaussModSolution{
    bool solvable;
    vector<int> x;
    vector<vector<int>> basis;
};

int expo(int a, int b, int mod){
    int res = 1;

    while (b){
        if (b & 1) res = (long long)res * a % mod;
        a = (long long)a * a % mod;
        b >>= 1;
    }

    return res;
}

/// target -= factor * source from column first on, (mod - factor) * source[l] + target[l] < 2^62 needs one % only
void subtract_row(vector<int>& target, const vector<int>& source, int factor, int first, int mod){
    long long negated = mod - factor;
    int width = target.size();
    for (int l = first; l < width; l++) target[l] = (target[l] + negated * source[l]) % mod;
}

/// Row echelon form of the first cols columns, pivots 1 in rows 0 to rank - 1, pos[j] = pivot row of column j or -1
vector<int> echelon_mod(vector<vector<int>>& a, int cols, int mod){
    int n = a.size(), i = 0;
    for (auto& row : a){
        for (auto& v : row){
            v %= mod;
            if (v < 0) v += mod;
        }
    }

    vector<int> pos(cols, -1);
    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        while (p < n && !a[p][j]) p++;
        if (p == n) continue;

        swap(a[p], a[i]);
        int width = a[i].size(), inv = expo(a[i][j], mod - 2, mod);
        for (int l = j; l < width; l++) a[i][l] = (long long)a[i][l] * inv % mod;

        for (int k = i + 1; k < n; k++){
            if (a[k][j]) subtract_row(a[k], a[i], a[k][j], j, mod);
        }

        pos[j] = i++;
    }

    return pos;
}

/// Reduced row echelon form: echelon form, then every pivot clears its column above it
vector<int> rref_mod(vector<vector<int>>& a, int cols, int mod){
    vector<int> pos = echelon_mod(a, cols, mod);
    for (int j = cols - 1; j >= 0; j--){
        if (pos[j] == -1) continue;

        int i = pos[j];
        for (int k = 0; k < i; k++){
            if (a[k][j]) subtract_row(a[k], a[i], a[k][j], j, mod);
        }
    }

    return pos;
}

GaussModSolution gauss(vector<vector<int>> equations, int mod){
    int n = equations.size(), m = equations[0].size() - 1;
    vector<int> pos = rref_mod(equations, m, mod);

    /// Free variable f at 1 and the other free ones at 0 force every pivot variable j to -equations[pos[j]][f]
    vector<vector<int>> basis;
    for (int f = 0; f < m; f++){
        if (pos[f] != -1) continue;

        vector<int> direction(m, 0);
        direction[f] = 1;
        for (int j = 0; j < f; j++){
            if (pos[j] != -1 && equations[pos[j]][f]) direction[j] = mod - equations[pos[j]][f];
        }
        basis.push_back(direction);
    }

    /// Only rows without a pivot can be inconsistent, as 0 = nonzero
    int rank = m - basis.size();
    for (int k = rank; k < n; k++){
        if (equations[k][m]) return {false, {}, basis};
    }

    vector<int> x(m, 0);
    for (int j = 0; j < m; j++){
        if (pos[j] != -1) x[j] = equations[pos[j]][m];
    }

    return {true, x, basis};
}

optional<vector<vector<int>>> matrix_inverse(const vector<vector<int>>& a, int mod){
    int n = a.size();
    vector<vector<int>> aug(n, vector<int>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = rref_mod(aug, n, mod);
    if (n && pos[n - 1] != n - 1) return nullopt;

    vector<vector<int>> inv(n, vector<int>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j];
    }

    return inv;
}

int matrix_rank(vector<vector<int>> a, int mod){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = echelon_mod(a, cols, mod);
    return cols - count(pos.begin(), pos.end(), -1);
}

int main(){
    /// Modulo 7, 2x + y - z = 1, -3x - y + 2z = 3 and -2x + y + 2z = 4 have the unique solution x = 2, y = 3, z = 6
    const int mod = 7;
    vector<vector<int>> equations = {{2, 1, -1, 1}, {-3, -1, 2, 3}, {-2, 1, 2, 4}};
    auto [solvable, x, basis] = gauss(equations, mod);
    assert(solvable);
    assert(basis.empty());
    assert((x == vector<int>{2, 3, 6}));

    /// x + 2y + 3z = 6 gives x = 6 - 2y - 3z, so the kernel is spanned by (-2, 1, 0) and (-3, 0, 1)
    GaussModSolution plane = gauss({{1, 2, 3, 6}}, mod);
    assert(plane.solvable);
    assert((plane.x == vector<int>{6, 0, 0}));
    assert((plane.basis == vector<vector<int>>{{5, 1, 0}, {4, 0, 1}}));

    /// x + y = 1 and x + y = 2 contradict each other, the kernel is still spanned by (-1, 1)
    GaussModSolution contradiction = gauss({{1, 1, 1}, {1, 1, 2}}, mod);
    assert(!contradiction.solvable);
    assert(contradiction.x.empty());
    assert((contradiction.basis == vector<vector<int>>{{6, 1}}));

    /// det = -2 is a unit modulo 7 but 0 modulo 2
    assert(matrix_rank({{1, 2}, {3, 4}}, mod) == 2);
    assert(matrix_rank({{1, 2}, {3, 4}}, 2) == 1);

    /// The inverse is -1/2 * {{4, -2}, {-3, 1}}, and 1/2 = 4 modulo 7
    optional<vector<vector<int>>> inv = matrix_inverse({{1, 2}, {3, 4}}, mod);
    assert(inv.has_value());
    assert((*inv == vector<vector<int>>{{5, 1}, {5, 3}}));

    assert(!matrix_inverse({{1, 2}, {3, 4}}, 2).has_value());

    return 0;
}
