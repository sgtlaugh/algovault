/***
 *
 * Gauss-Jordan Elimination modulo a prime
 * Linear systems, solution spaces, rank and inverse of matrices modulo a prime
 *
 * Complexity: O(min(n, m) * n * m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * mod must be a prime below 2^31, entries may be any int and are reduced into [0, mod)
 *
 * gauss(equations, res, mod): n >= 1 equations in m variables, equations[i][m] is the right-hand side
 *     returns -1 when there is no solution
 *     otherwise returns the number of free variables f, 0 for a unique solution, and res holds one solution
 * gauss(equations, res, basis, mod): as above, and basis holds f vectors spanning the kernel {x : Ax = 0}
 *     every solution is res + c_1 * basis[0] + ... + c_f * basis[f - 1] for exactly one choice of the c_k
 *     basis is filled even when there is no solution, it depends only on the coefficients
 * matrix_rank(a, mod): rank of any n x m matrix
 * matrix_inverse(a, inv, mod): for a square a, returns false when it is singular, otherwise inv = a^-1
 *
 * For instance, the system of linear equations modulo 7
 *
 *     2x + y - z   = 1  ----->  equations[0] = {2, 1, 6, 1}
 *     -3x - y + 2z = 3  ----->  equations[1] = {4, 6, 2, 3}
 *     -2x + y + 2z = 4  ----->  equations[2] = {5, 1, 2, 4}
 *
 * has the unique solution x = 2, y = 3, z = 6
 *
 * For problems on graphs, make sure the graph is connected and a single component
 * If not, then re-number the vertices and solve for each component separately
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

int expo(int a, int b, int mod){
    int res = 1;

    while (b){
        if (b & 1) res = (long long)res * a % mod;
        a = (long long)a * a % mod;
        b >>= 1;
    }

    return res;
}

/// Brings a into reduced row echelon form over its first cols columns, pivots become 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
vector<int> rref_mod(vector<vector<int>>& a, int cols, int mod){
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

        for (int k = 0; k < n; k++){
            if (k == i || !a[k][j]) continue;
            int x = a[k][j];
            for (int l = j; l < width; l++){
                if (!a[i][l]) continue;
                long long v = (a[k][l] - (long long)a[i][l] * x) % mod;
                a[k][l] = v < 0 ? v + mod : v;
            }
        }

        pos[j] = i++;
    }

    return pos;
}

int gauss(vector<vector<int>> equations, vector<int>& res, vector<vector<int>>& basis, int mod){
    int n = equations.size(), m = equations[0].size() - 1;
    vector<int> pos = rref_mod(equations, m, mod);

    basis.clear();
    for (int f = 0; f < m; f++){
        if (pos[f] != -1) continue;
        vector<int> x(m, 0);
        x[f] = 1;
        for (int j = 0; j < f; j++){
            if (pos[j] != -1 && equations[pos[j]][f]) x[j] = mod - equations[pos[j]][f];
        }
        basis.push_back(x);
    }

    res.assign(m, 0);
    int rank = m - basis.size();
    for (int k = rank; k < n; k++){
        if (equations[k][m]) return -1;
    }

    for (int j = 0; j < m; j++){
        if (pos[j] != -1) res[j] = equations[pos[j]][m];
    }

    return basis.size();
}

int gauss(vector<vector<int>> equations, vector<int>& res, int mod){
    vector<vector<int>> basis;
    return gauss(equations, res, basis, mod);
}

bool matrix_inverse(const vector<vector<int>>& a, vector<vector<int>>& inv, int mod){
    int n = a.size();
    vector<vector<int>> aug(n, vector<int>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = rref_mod(aug, n, mod);
    if (n && pos[n - 1] != n - 1) return false;

    inv.assign(n, vector<int>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j];
    }

    return true;
}

int matrix_rank(vector<vector<int>> a, int mod){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = rref_mod(a, cols, mod);
    return cols - count(pos.begin(), pos.end(), -1);
}

const int MOD = 1000000007;

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;

    for (int tc = 1; tc <= t; tc++){
        int n, m, k;
        if (scanf("%d %d %d", &n, &m, &k) != 3) return 0;

        /// row v encodes sum of neighbour colours - color(v) = 0 (mod K), every free variable takes any of K colours
        vector<vector<int>> a(n, vector<int>(n, 0));
        for (int v = 0; v < n; v++) a[v][v] = -1;
        for (int i = 0; i < m; i++){
            int u, v;
            if (scanf("%d %d", &u, &v) != 2) return 0;
            u--, v--;
            a[u][v] = a[v][u] = 1;
        }

        int free_vars = n - matrix_rank(a, k);
        printf("Case %d: %d\n", tc, expo(k, free_vars, MOD));
    }
    return 0;
}
