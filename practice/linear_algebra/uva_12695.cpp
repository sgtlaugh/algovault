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

            vector<double> p = gauss(eq).x;
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
