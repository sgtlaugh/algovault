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

const int dx[] = {-1, 1, 0, 0};
const int dy[] = {0, 0, -1, 1};

int w, h;
vector<string> grid;
vector<vector<bool>> reach;

bool inside(int x, int y){
    return x >= 0 && x < h && y >= 0 && y < w && grid[x][y] != 'X';
}

bool live(int x, int y){
    return x >= 0 && x < h && grid[x][y] == '.' && reach[x][y];
}

int main(){
    int u, d, l, r;

    for (int tc = 0; scanf("%d %d", &w, &h) == 2 && (w || h); tc++){
        scanf("%d %d %d %d", &u, &d, &l, &r);
        double p[4] = {u / 100.0, d / 100.0, l / 100.0, r / 100.0};
        grid.assign(h, "");
        char buf[32];
        for (int i = 0; i < h; i++){
            scanf("%s", buf);
            grid[i] = buf;
        }

        reach.assign(h, vector<bool>(w, false));
        vector<pair<int, int>> stk;
        int tops = 0;
        for (int j = 0; j < w; j++){
            if (grid[0][j] == '.') reach[0][j] = true, stk.push_back({0, j}), tops++;
        }
        while (!stk.empty()){
            int x = stk.back().first, y = stk.back().second;
            stk.pop_back();
            for (int k = 0; k < 4; k++){
                int nx = x + dx[k], ny = y + dy[k];
                if (p[k] > 0 && inside(nx, ny) && grid[nx][ny] == '.' && !reach[nx][ny]){
                    reach[nx][ny] = true;
                    stk.push_back({nx, ny});
                }
            }
        }

        vector<vector<vector<double>>> sinv(h);
        vector<vector<double>> g(h, vector<double>(w, 0)), up(h, vector<double>(w, 0));
        for (int x = 0; x < h; x++){
            vector<vector<double>> a(w, vector<double>(w, 0));
            vector<double> down(w, 0);
            for (int y = 0; y < w; y++){
                a[y][y] += 1;
                if (!live(x, y)) continue;
                if (x == 0) g[x][y] = 1.0 / tops;
                if (live(x - 1, y)) down[y] = p[1];
                if (live(x + 1, y)) up[x][y] = p[0];
                for (int k = 0; k < 4; k++){
                    int nx = x + dx[k], ny = y + dy[k];
                    if (!inside(nx, ny)) a[y][y] -= p[k];
                    else if (nx == x && grid[nx][ny] == '.') a[ny][y] -= p[k];
                }
            }

            if (x > 0){
                for (int i = 0; i < w; i++){
                    if (down[i] == 0) continue;
                    for (int j = 0; j < w; j++){
                        a[i][j] -= down[i] * sinv[x - 1][i][j] * up[x - 1][j];
                        g[x][i] += down[i] * sinv[x - 1][i][j] * g[x - 1][j];
                    }
                }
            }

            matrix_inverse(a, sinv[x], 1e-14);
        }

        vector<vector<double>> v(h, vector<double>(w, 0));
        for (int x = h - 1; x >= 0; x--){
            vector<double> rhs = g[x];
            if (x + 1 < h){
                for (int j = 0; j < w; j++) rhs[j] += up[x][j] * v[x + 1][j];
            }
            for (int i = 0; i < w; i++){
                for (int j = 0; j < w; j++) v[x][i] += sinv[x][i][j] * rhs[j];
            }
        }

        if (tc) puts("");
        for (int x = 0; x < h; x++){
            for (int y = 0; y < w; y++){
                if (grid[x][y] != 'T') continue;
                double prob = 0;
                for (int k = 0; k < 4; k++){
                    int px = x - dx[k], py = y - dy[k];
                    if (py >= 0 && py < w && live(px, py)) prob += v[px][py] * p[k];
                }
                printf("%.9f\n", prob);
            }
        }
    }
    return 0;
}
