/***
 *
 * Gauss-Jordan Elimination
 * Linear systems, rank and inverse of matrices over the reals
 *
 * Complexity: O(min(n, m) n m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, eps = 1e-9): {solvable, x, free_variables} for n >= 1 equations in m variables, equations[i][m] is the right-hand side
 * matrix_rank(a, eps = 1e-9): rank of any n x m matrix
 * matrix_inverse(a, eps = 1e-9): a^-1 for a square a, nullopt when a is singular
 *
 * x is one solution, empty when solvable is false, and free_variables is m minus the rank either way
 * eps is relative, scaled by max(1, largest |coefficient| of the input), a pivot within it counts as 0
 * The right-hand side of gauss is left out of that scale, so x = 1e9 with coefficients of 1 still solves
 * A row of gauss that reduces to 0 = r is consistent when |r| is within eps scaled by the largest |entry| including it
 * T may be an exact type such as a fraction with eps = 0, and long double or __float128 give more precision
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <class T>
struct GaussSolution{
    bool solvable;
    vector<T> x;
    int free_variables;
};

/// Relative, an absolute eps rejects solvable systems with large coefficients
/// Over the first cols columns, a large right-hand side in the scale would push ordinary pivots under it
template <class T>
T tolerance(const vector<vector<T>>& a, int cols, T eps){
    T tol = 1;
    for (auto& row : a) for (int j = 0; j < cols; j++) tol = max(tol, (T)abs(row[j]));
    return tol * eps;
}

/// Brings the first cols columns of a to row echelon form with partial pivoting, pivots are not scaled to 1
/// Returns pos, pos[j] = row of the pivot in column j or -1, pivot rows are 0, 1, ..., rank - 1 in column order
/// Only rows below a pivot are cleared, about 2/3 of the work of clearing above too
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
        for (int k = i + 1; k < n; k++){
            if (a[k][j] == 0) continue;
            T factor = a[k][j] / a[i][j];
            for (int l = j; l < width; l++) a[k][l] -= a[i][l] * factor;
        }

        pos[j] = i++;
    }

    return pos;
}

template <class T>
GaussSolution<T> gauss(vector<vector<T>> equations, const T eps=1e-9){
    int n = equations.size(), m = equations[0].size() - 1;
    T tol = tolerance(equations, m + 1, eps);
    vector<int> pos = row_reduce(equations, m, tolerance(equations, m, eps));
    int rank = m - count(pos.begin(), pos.end(), -1);

    /// Only rows without a pivot can be inconsistent, as 0 = nonzero
    for (int k = rank; k < n; k++){
        if (abs(equations[k][m]) > tol) return {false, {}, m - rank};
    }

    /// Back substitution from the last pivot up, free variables stay 0
    vector<T> x(m, 0);
    for (int j = m - 1; j >= 0; j--){
        if (pos[j] == -1) continue;

        vector<T>& row = equations[pos[j]];
        T value = row[m];
        for (int l = j + 1; l < m; l++) value -= row[l] * x[l];
        x[j] = value / row[j];
    }

    return {true, x, m - rank};
}

template <class T>
int matrix_rank(vector<vector<T>> a, const T eps=1e-9){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = row_reduce(a, cols, tolerance(a, cols, eps));
    return cols - count(pos.begin(), pos.end(), -1);
}

template <class T>
optional<vector<vector<T>>> matrix_inverse(const vector<vector<T>>& a, const T eps=1e-9){
    int n = a.size();
    vector<vector<T>> aug(n, vector<T>(2 * n, 0));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) aug[i][j] = a[i][j];
        aug[i][n + i] = 1;
    }

    vector<int> pos = row_reduce(aug, n, tolerance(a, n, eps));
    if (n && pos[n - 1] != n - 1) return nullopt;

    /// Clears above each pivot from the last up, the left half of a pivot row is then 0 past its pivot
    for (int i = n - 1; i >= 0; i--){
        for (int k = 0; k < i; k++){
            if (aug[k][i] == 0) continue;
            T factor = aug[k][i] / aug[i][i];
            for (int l = n; l < 2 * n; l++) aug[k][l] -= aug[i][l] * factor;
        }
    }

    vector<vector<T>> inv(n, vector<T>(n));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) inv[i][j] = aug[i][n + j] / aug[i][i];
    }

    return inv;
}

int main(){
    int n, m, a, b;
    scanf("%d %d %d %d", &n, &m, &a, &b);
    vector<vector<int>> adj(n);
    for (int i = 0; i < m; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u - 1].push_back(v - 1);
        adj[v - 1].push_back(u - 1);
    }

    vector<double> stay(n);
    for (auto& p : stay) scanf("%lf", &p);

    /// moves[i] lists every room reachable from i in one minute with its probability
    vector<vector<pair<int, double>>> moves(n);
    for (int i = 0; i < n; i++){
        moves[i].push_back({i, stay[i]});
        for (int j : adj[i]) moves[i].push_back({j, (1 - stay[i]) / adj[i].size()});
    }

    /// visits[(i, j)] = expected visits to the pair, which for i == j is the chance they meet in room i
    int states = n * n, start = (a - 1) * n + (b - 1);
    vector<vector<double>> equations(states, vector<double>(states + 1, 0));
    for (int s = 0; s < states; s++) equations[s][s] = 1;
    equations[start][states] = 1;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j) continue;
            for (auto [ni, pi] : moves[i]){
                for (auto [nj, pj] : moves[j]) equations[ni * n + nj][i * n + j] -= pi * pj;
            }
        }
    }

    vector<double> visits = gauss(equations).x;
    for (int i = 0; i < n; i++) printf("%.10f ", visits[i * n + i]);
    puts("");
    return 0;
}
