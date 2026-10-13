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

int main(){
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, m;
        scanf("%d %d", &n, &m);
        vector<int> out_deg(n, 0);
        vector<vector<vector<int>>> cnt(26, vector<vector<int>>(n, vector<int>(n, 0)));
        vector<vector<bool>> road(n, vector<bool>(n, false));
        while (m--){
            int u, v;
            char w[4];
            scanf("%d %d %s", &u, &v, w);
            cnt[w[0] - 'A'][u][v]++, out_deg[u]++, road[u][v] = true;
        }

        /// p(a, b): total chance of the walks from a to b that read a palindrome, pairing the first road out of a
        /// with the last road into b. step[a][b][x][y] weighs pairing a -> x with y -> b, direct[a][b] is a single
        /// middle road a -> b
        auto prob = [&](int u, int v, int c){ return (double)cnt[c][u][v] / out_deg[u]; };
        vector<vector<double>> direct(n, vector<double>(n, 0));
        vector<vector<vector<vector<double>>>> step(n, vector<vector<vector<double>>>(n, vector<vector<double>>(n, vector<double>(n, 0))));
        for (int a = 0; a < n; a++){
            for (int b = 0; b < n; b++){
                for (int c = 0; c < 26; c++){
                    direct[a][b] += prob(a, b, c);
                    for (int x = 0; x < n; x++){
                        if (!cnt[c][a][x]) continue;
                        for (int y = 0; y < n; y++){
                            if (cnt[c][y][b]) step[a][b][x][y] += prob(a, x, c) * prob(y, b, c);
                        }
                    }
                }
            }
        }

        /// The walk stops at t and in cities that cannot reach t, so only the other cities that reach t may sit
        /// strictly inside a path; keeping the rest out also keeps the system nonsingular
        vector<vector<double>> answer(n, vector<double>(n, 0));
        for (int e = 0; e < n; e++){
            vector<bool> inside(n, false);
            vector<int> queue = {e};
            for (int i = 0; i < (int)queue.size(); i++){
                for (int u = 0; u < n; u++){
                    if (road[u][queue[i]] && !inside[u] && u != e) inside[u] = true, queue.push_back(u);
                }
            }

            vector<int> cities, pos(n, -1);
            for (int u = 0; u < n; u++){
                if (inside[u]) pos[u] = cities.size(), cities.push_back(u);
            }
            int k = cities.size(), vars = k * k;
            if (!k) continue;

            /// p(a, a) also counts the empty walk: meeting in a city may close an even palindrome or go on
            /// with a palindromic loop back to it
            auto expand = [&](int a, int b, vector<double>& row){
                for (int x : cities){
                    for (int y : cities) row[pos[x] * k + pos[y]] -= step[a][b][x][y];
                }
                return direct[a][b] + (a == b);
            };

            vector<vector<double>> eq(vars, vector<double>(vars + 1, 0));
            for (int a : cities){
                for (int b : cities){
                    vector<double>& row = eq[pos[a] * k + pos[b]];
                    row[pos[a] * k + pos[b]] += 1;
                    row[vars] = expand(a, b, row);
                }
            }

            vector<double> p;
            gauss(eq, p);
            for (int s : cities){
                vector<double> row(vars, 0);
                double res = expand(s, e, row);
                for (int i = 0; i < vars; i++) res -= row[i] * p[i];
                answer[s][e] = res;
            }
        }

        printf("Case %d:\n", tc);
        int q;
        scanf("%d", &q);
        while (q--){
            int s, e;
            scanf("%d %d", &s, &e);
            printf("%.6f\n", answer[s][e] + 1e-9);
        }
    }
    return 0;
}
