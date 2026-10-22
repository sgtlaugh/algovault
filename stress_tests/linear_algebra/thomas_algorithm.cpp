// LINK: -lquadmath
#include "../common.h"
#include <quadmath.h>

#define main library_main
#include "../../code_library/linear_algebra/thomas_algorithm.cpp"
#undef main

/// Dense elimination with partial pivoting in quad precision, sharing nothing with the tridiagonal sweep
vector<__float128> dense_solve(vector<vector<__float128>> a, vector<__float128> b){
    int n = b.size();
    for (int c = 0; c < n; c++){
        int p = c;
        for (int r = c + 1; r < n; r++) if (fabsq(a[r][c]) > fabsq(a[p][c])) p = r;
        swap(a[p], a[c]), swap(b[p], b[c]);
        for (int r = c + 1; r < n; r++){
            __float128 f = a[r][c] / a[c][c];
            for (int k = c; k < n; k++) a[r][k] -= f * a[c][k];
            b[r] -= f * b[c];
        }
    }

    vector<__float128> x(n);
    for (int r = n - 1; r >= 0; r--){
        __float128 s = b[r];
        for (int k = r + 1; k < n; k++) s -= a[r][k] * x[k];
        x[r] = s / a[r][r];
    }

    return x;
}

int main(){
    for (long long it = 0; it < stress::scaled(8000); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 3000);
        bool weak = it % 3 == 0;  /// the header's weak dominance case: rows with equality, a strict first row, every lower[i] nonzero past the first
        vector<long double> lower(n), diagonal(n), upper(n), rhs(n);
        for (int i = 0; i < n; i++){
            /// Coefficients outside the matrix must be ignored, so they get junk
            lower[i] = stress::rand_int(-1000, 1000) / 100.0L, upper[i] = stress::rand_int(-1000, 1000) / 100.0L;
            if (weak && i && lower[i] == 0) lower[i] = 1;

            long double off = fabsl(i ? lower[i] : 0) + fabsl(i + 1 < n ? upper[i] : 0);
            long double extra = weak && i ? (stress::rand_int(0, 1) ? 0 : stress::rand_int(1, 1000) / 100.0L) : stress::rand_int(1, 1000) / 100.0L;
            diagonal[i] = (off + extra) * (stress::rand_int(0, 1) ? 1 : -1);
            rhs[i] = stress::rand_int(-100000, 100000) / 64.0L;
        }

        auto x = solve_tridiagonal(lower, diagonal, upper, rhs);
        assert((int)x.size() == n);

        long double norm = 1;
        for (auto v : x) norm = max(norm, fabsl(v));
        for (int i = 0; i < n; i++){
            long double s = diagonal[i] * x[i] + (i ? lower[i] * x[i - 1] : 0) + (i + 1 < n ? upper[i] * x[i + 1] : 0);
            assert(fabsl(s - rhs[i]) <= 1e-17L * norm * (fabsl(diagonal[i]) + fabsl(lower[i]) + fabsl(upper[i])));
        }

        /// Weak dominance can be badly conditioned, only its residual above is bounded
        if (n <= 60 && !weak){
            vector<vector<__float128>> a(n, vector<__float128>(n, 0));
            vector<__float128> b(n);
            for (int i = 0; i < n; i++){
                a[i][i] = diagonal[i];
                if (i) a[i][i - 1] = lower[i];
                if (i + 1 < n) a[i][i + 1] = upper[i];
                b[i] = rhs[i];
            }
            auto expected = dense_solve(a, b);
            for (int i = 0; i < n; i++) assert(fabsq(expected[i] - x[i]) <= (__float128)1e-17L * norm);
        }
    }

    /// M-matrix version on well-conditioned random systems: an excess of at least 1/2 in every row, checked per component against quad
    for (long long it = 0; it < stress::scaled(8000); it++){
        int n = stress::rand_int(0, it % 10 ? 12 : 60);
        vector<long double> lower(n), upper(n), excess(n), rhs(n);
        for (int i = 0; i < n; i++){
            lower[i] = stress::rand_int(0, 3) ? stress::rand_int(0, 1000) / 100.0L : 0;
            upper[i] = stress::rand_int(0, 3) ? stress::rand_int(0, 1000) / 100.0L : 0;
            excess[i] = stress::rand_int(50, 500) / 100.0L;
            rhs[i] = stress::rand_int(0, 3) ? stress::rand_int(0, 100000) / 64.0L : 0;
        }

        auto x = solve_tridiagonal_m_matrix(lower, upper, excess, rhs);
        assert((int)x.size() == n);
        if (n == 0) continue;

        vector<vector<__float128>> a(n, vector<__float128>(n, 0));
        vector<__float128> b(n);
        for (int i = 0; i < n; i++){
            a[i][i] = (__float128)lower[i] + upper[i] + excess[i];  /// lower[0] and upper[n - 1] count in the diagonal
            if (i) a[i][i - 1] = -lower[i];
            if (i + 1 < n) a[i][i + 1] = -upper[i];
            b[i] = rhs[i];
        }
        auto expected = dense_solve(a, b);
        __float128 largest = 0;
        for (auto v : expected) largest = max(largest, fabsq(v));
        /// Pivoted quad elimination is only accurate relative to the largest entry, so the floor scales with it
        for (int i = 0; i < n; i++) assert(fabsq(expected[i] - x[i]) <= (__float128)1e-16L * fabsq(expected[i]) + (__float128)1e-30L * largest);
    }

    /// Walks with strong drift, where elimination with subtraction loses every digit: cell i steps right with chance p[i]
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 10 ? 50 : 400);
        vector<long double> p(n);
        int style = stress::rand_int(0, 2);
        for (int i = 0; i < n; i++){
            if (style == 0) p[i] = stress::rand_int(1, 99) / 100.0L;
            else if (style == 1) p[i] = 2 * i < n ? 0.99L : 0.01L;  /// a trap in the middle, as in UVa 1413
            else p[i] = stress::rand_int(0, 1) ? 0.99L : stress::rand_int(1, 99) / 100.0L;
        }

        vector<long double> lower(n), upper(n), excess(n, 0), rhs(n, 0);
        for (int i = 0; i < n; i++) lower[i] = 1 - p[i], upper[i] = p[i];
        rhs[n - 1] = p[n - 1];
        auto h = solve_tridiagonal_m_matrix(lower, upper, excess, rhs);

        /// The chance to leave on the right first is a ratio of sums of positive products, d_0 = 1 and d_j = d_(j - 1) (1 - p[j]) / p[j]
        vector<long double> prefix(n + 1, 0);
        long double d = 1;
        for (int j = 0; j <= n; j++){
            if (j) d *= (1 - p[j - 1]) / p[j - 1];
            prefix[j] = (j ? prefix[j - 1] : 0) + d;
        }
        for (int i = 0; i < n; i++){
            long double expected = prefix[i] / prefix[n];
            assert(fabsl(h[i] - expected) <= 1e-15L * expected);
        }
    }

    return 0;
}
