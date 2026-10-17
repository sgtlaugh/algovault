/***
 *
 * Gauss-Jordan Elimination modulo a prime
 * Linear systems, solution spaces, rank and inverse of matrices modulo a prime
 *
 * Complexity: O(min(n, m) n m + m^2) for gauss, O(min(n, m) n m) for matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, mod): {solvable, x, basis} for n >= 1 equations in m variables, equations[i][m] is the right-hand side
 * matrix_rank(a, mod): rank of any n x m matrix
 * matrix_inverse(a, mod): a^-1 for a square a, nullopt when a is singular
 *
 * x is one solution, empty when solvable is false
 * basis spans the kernel {x : Ax = 0}, every solution is x + c_1 basis[0] + ... + c_f basis[f - 1] for exactly one choice of the c_k
 * basis is filled even when there is no solution, it depends only on the coefficients
 * mod must be a prime below 2^31, entries may be any int and are reduced into [0, mod)
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

/// target -= factor * source over the columns from first on
/// (mod - factor) * source[l] + target[l] stays below 2^62 for mod < 2^31, so one % and no sign fix
void subtract_row(vector<int>& target, const vector<int>& source, int factor, int first, int mod){
    long long negated = mod - factor;
    int width = target.size();
    for (int l = first; l < width; l++) target[l] = (target[l] + negated * source[l]) % mod;
}

/// Brings a into row echelon form over its first cols columns, pivots become 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
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

int matrix_rank(vector<vector<int>> a, int mod){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = echelon_mod(a, cols, mod);
    return cols - count(pos.begin(), pos.end(), -1);
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

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<int>> equations(n, vector<int>(m + 1));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++) scanf("%d", &equations[i][j]);
    }
    for (int i = 0; i < n; i++) scanf("%d", &equations[i][m]);

    auto [solvable, x, basis] = gauss(equations, 998244353);
    if (!solvable){
        puts("-1");
        return 0;
    }
    printf("%d\n", (int)basis.size());
    for (int j = 0; j < m; j++) printf("%d%c", x[j], j + 1 == m ? '\n' : ' ');
    for (auto& direction : basis){
        for (int j = 0; j < m; j++) printf("%d%c", direction[j], j + 1 == m ? '\n' : ' ');
    }
    return 0;
}
