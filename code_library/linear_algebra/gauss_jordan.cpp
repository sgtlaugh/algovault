/***
 *
 * Gauss-Jordan Elimination
 * Linear systems, rank and inverse of matrices over the reals
 *
 * Complexity: O(min(n, m) n m) for gauss and matrix_rank, O(n^3) for matrix_inverse
 *
 * gauss(equations, eps = 1e-9): {solvable, x, free_variables} for n >= 1 equations in m variables, the right side last
 * matrix_inverse(a, eps = 1e-9): a^-1 for a square a, nullopt when a is singular
 * matrix_rank(a, eps = 1e-9): rank of any n x m matrix
 *
 * Notes:
 *   - x is one solution, empty when solvable is false, and free_variables is m minus the rank either way
 *   - eps is relative, scaled by max(1, largest |coefficient| of the input), a pivot within it counts as 0
 *   - The right-hand side of gauss is left out of that scale, so x = 1e9 with coefficients of 1 still solves
 *   - A gauss row reducing to 0 = r is consistent when |r| is within eps scaled by the largest |entry| including it
 *   - T may be an exact type such as a fraction with eps = 0, and long double or __float128 give more precision
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

/// Relative, as an absolute eps rejects large coefficients, and over cols columns, as a large right side hides pivots
template <class T>
T tolerance(const vector<vector<T>>& a, int cols, T eps){
    T tol = 1;
    for (auto& row : a) for (int j = 0; j < cols; j++) tol = max(tol, (T)abs(row[j]));
    return tol * eps;
}

/// Row echelon form of the first cols columns, unscaled pivots in rows 0 to rank - 1, pos[j] = pivot row of j or -1
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

template <class T>
int matrix_rank(vector<vector<T>> a, const T eps=1e-9){
    int cols = a.empty() ? 0 : a[0].size();
    vector<int> pos = row_reduce(a, cols, tolerance(a, cols, eps));
    return cols - count(pos.begin(), pos.end(), -1);
}

int main(){
    /// 2x + y - z = 8, -3x - y + 2z = -11 and -2x + y + 2z = -3 have the unique solution x = 2, y = 3, z = -1
    vector<vector<double>> equations = {{2, 1, -1, 8}, {-3, -1, 2, -11}, {-2, 1, 2, -3}};
    auto [solvable, x, free_variables] = gauss(equations);
    assert(solvable);
    assert(free_variables == 0);
    assert(abs(x[0] - 2) < 1e-9);
    assert(abs(x[1] - 3) < 1e-9);
    assert(abs(x[2] + 1) < 1e-9);

    /// Without the third equation z is free and stays 0, then 2x + y = 8 and -3x - y = -11 give x = 3, y = 2
    GaussSolution<double> partial = gauss(vector<vector<double>>{{2, 1, -1, 8}, {-3, -1, 2, -11}});
    assert(partial.solvable);
    assert(partial.free_variables == 1);
    assert(abs(partial.x[0] - 3) < 1e-9);
    assert(abs(partial.x[1] - 2) < 1e-9);
    assert(abs(partial.x[2]) < 1e-9);

    /// x + y = 2 and 2x + 2y = 5 contradict each other
    GaussSolution<double> contradiction = gauss(vector<vector<double>>{{1, 1, 2}, {2, 2, 5}});
    assert(!contradiction.solvable);
    assert(contradiction.x.empty());

    /// The third row is twice the second minus the first
    assert(matrix_rank(vector<vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}) == 2);

    /// det = -2, so the inverse is -1/2 * {{4, -2}, {-3, 1}}
    optional<vector<vector<double>>> inv = matrix_inverse(vector<vector<double>>{{1, 2}, {3, 4}});
    assert(inv.has_value());
    assert(abs((*inv)[0][0] + 2) < 1e-9);
    assert(abs((*inv)[0][1] - 1) < 1e-9);
    assert(abs((*inv)[1][0] - 1.5) < 1e-9);
    assert(abs((*inv)[1][1] + 0.5) < 1e-9);

    assert(!matrix_inverse(vector<vector<double>>{{1, 2}, {2, 4}}).has_value());

    return 0;
}
