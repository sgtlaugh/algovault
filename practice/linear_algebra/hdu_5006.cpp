/***
 *
 * Gauss-Jordan Elimination
 * Linear systems, rank and inverse of matrices over the reals
 *
 * Complexity: O(min(n, m) * n * m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, res): n >= 1 equations in m variables, equations[i][m] is the right-hand side
 *     returns -1 when there is no solution
 *     otherwise returns the number of free variables, 0 for a unique solution, and res holds one solution
 * matrix_rank(a): rank of any n x m matrix
 * matrix_inverse(a, inv): for a square a, returns false when it is singular, otherwise inv = a^-1
 *
 * For instance, the system of linear equations
 *
 *     2x + y - z   = 8    ----->  equations[0] = {2, 1, -1, 8}
 *     -3x - y + 2z = -11  ----->  equations[1] = {-3, -1, 2, -11}
 *     -2x + y + 2z = -3   ----->  equations[2] = {-2, 1, 2, -3}
 *
 * has the unique solution x = 2, y = 3, z = -1
 *
 * For problems on graphs, make sure the graph is connected and a single component
 * If not, then re-number the vertices and solve for each component separately
 *
 * Every function takes an optional eps, scaled by max(1, largest |coefficient| of its input)
 * A pivot within that tolerance of 0 counts as 0, making the matrix rank deficient or singular
 * The right-hand side of gauss is left out of the pivot scale, so x = 1e9 with coefficients of 1 still solves
 * A row of gauss that reduces to 0 = r is consistent when |r| is within eps scaled by the largest |entry| including it
 * If more precision is required, use long double or __float128 from quadmath.h if supported
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Relative, an absolute eps rejects solvable systems with large coefficients
/// Over the first cols columns, a large right-hand side in the scale would push ordinary pivots under it
template <class T>
T tolerance(const vector<vector<T>>& a, int cols, T eps){
    T tol = 1;
    for (auto& row : a) for (int j = 0; j < cols; j++) tol = max(tol, (T)abs(row[j]));
    return tol * eps;
}

/// Eliminates the first cols columns of a with partial pivoting, pivots are not scaled to 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
template <class T>
vector<int> row_reduce(vector<vector<T>>& a, int cols, T tol){
    int n = a.size(), i = 0;
    vector<int> pos(cols, -1);

    for (int j = 0; j < cols && i < n; j++){
        int p = i;
        for (int k = i; k < n; k++){
            if (abs(a[k][j]) > abs(a[p][j])) p = k;
        }
        if (abs(a[p][j]) <= tol) continue;

        swap(a[p], a[i]);
        int width = a[i].size();
        for (int k = 0; k < n; k++){
            if (k == i || a[k][j] == 0) continue;
            T x = a[k][j] / a[i][j];
            for (int l = j; l < width; l++) a[k][l] -= a[i][l] * x;
        }

        pos[j] = i++;
    }

    return pos;
}

template <class T>
int gauss(vector<vector<T>> equations, vector<T>& res, const T eps=1e-9){
    int n = equations.size(), m = equations[0].size() - 1, f_var = 0;
    T tol = tolerance(equations, m + 1, eps);
    vector<int> pos = row_reduce(equations, m, tolerance(equations, m, eps));

    res.assign(m, 0);
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = nonzero
    for (int k = rank; k < n; k++){
        if (abs(equations[k][m]) > tol) return -1;
    }

    for (int j = 0; j < m; j++){
        if (pos[j] == -1) f_var++;
        else res[j] = equations[pos[j]][m] / equations[pos[j]][j];
    }

    return f_var;
}

template <class T>
bool matrix_inverse(const vector<vector<T>>& a, vector<vector<T>>& inv, const T eps=1e-9){
    int n = a.size();
    vector<vector<T>> aug(n, vector<T>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = row_reduce(aug, n, tolerance(a, n, eps));
    if (n && pos[n - 1] != n - 1) return false;

    inv.assign(n, vector<T>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j] / aug[i][i];
    }

    return true;
}

template <class T>
int matrix_rank(vector<vector<T>> a, const T eps=1e-9){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = row_reduce(a, cols, tolerance(a, cols, eps));
    return cols - count(pos.begin(), pos.end(), -1);
}

struct DSU{
    vector<int> par;

    DSU(int n) : par(n){
        iota(par.begin(), par.end(), 0);
    }

    int find(int x){
        while (par[x] != x) x = par[x] = par[par[x]];
        return x;
    }

    void unite(int a, int b){
        par[find(a)] = find(b);
    }
};

int main(){
    int t;
    scanf("%d", &t);
    while (t--){
        int n, m, s, e;
        scanf("%d %d %d %d", &n, &m, &s, &e);
        s--, e--;

        /// 0 ohm wires make their endpoints a single node, which shrinks the random graph to a few hundred nodes
        DSU zero(n), all(n);
        vector<pair<int, int>> ones;
        for (int i = 0; i < m; i++){
            int u, v, c;
            scanf("%d %d %d", &u, &v, &c);
            u--, v--;
            all.unite(u, v);
            if (c == 0) zero.unite(u, v);
            else ones.push_back({u, v});
        }

        if (all.find(s) != all.find(e)){
            puts("inf");
            continue;
        }
        if (zero.find(s) == zero.find(e)){
            puts("0.000000");
            continue;
        }

        /// Only the component holding s and t, any other component would leave the system singular
        vector<int> id(n, -1);
        int k = 0, src = zero.find(s), sink = zero.find(e);
        for (int i = 0; i < n; i++){
            int r = zero.find(i);
            if (all.find(i) == all.find(s) && r != sink && id[r] == -1) id[r] = k++;
        }

        /// 1 A enters at s and the potential of t is 0, so the potential of s is the resistance
        vector<vector<double>> equations(k, vector<double>(k + 1, 0));
        for (auto& w : ones){
            int a = zero.find(w.first), b = zero.find(w.second);
            if (a == b || all.find(a) != all.find(s)) continue;
            int x = id[a], y = id[b];
            if (x != -1){
                equations[x][x] += 1;
                if (y != -1) equations[x][y] -= 1;
            }
            if (y != -1){
                equations[y][y] += 1;
                if (x != -1) equations[y][x] -= 1;
            }
        }
        equations[id[src]][k] = 1;

        vector<double> res;
        gauss(equations, res);
        printf("%.6f\n", res[id[src]]);
    }
    return 0;
}
