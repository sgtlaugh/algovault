// LINK: -lquadmath
#include "common.h"
#include <quadmath.h>

#define main library_main
#include "../code_library/gauss_band_matrix.cpp"
#undef main

typedef GaussBand<long double> Band;

/// Dense elimination with partial pivoting in quad precision
vector<__float128> dense_solve(vector<vector<__float128>> a, vector<__float128> b){
    int n = b.size();
    for (int c = 0; c < n; c++){
        int p = c;
        for (int r = c + 1; r < n; r++) if (fabsq(a[r][c]) > fabsq(a[p][c])) p = r;
        swap(a[p], a[c]), swap(b[p], b[c]);
        for (int r = c + 1; r < n; r++){
            __float128 f = a[r][c] / a[c][c];
            if (f == 0) continue;
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
    /// Lazily committed pages, the struct reserves room for the largest grid
    Band* gauss = (Band*)operator new(sizeof(Band));

    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, it % 10 ? 12 : 400), m = stress::rand_int(1, it % 10 ? 8 : 6);
        if (it % 50 == 0) m = 31, n = min(n, 40);  /// the widest band, kept small: cells and band size multiply
        int band_size = it % 3 ? 2 * m + 3 : 2 * stress::rand_int(0, min(34, 2 * m + 2)) + 1, h = band_size / 2, cells = n * m;
        long double rhs_default = stress::rand_int(-5, 5);
        new (gauss) Band(n, m, it % 3 ? 0 : band_size, rhs_default);

        /// Cell numbers as the library assigns them, so that coefficients land inside the band
        auto cell = [&](int u){ return make_pair((cells - 1 - u) / m, (cells - 1 - u) % m); };

        vector<map<int, long double>> coef(cells);
        vector<long double> rhs(cells, rhs_default);
        for (int u = 0; u < cells; u++){
            long double off = 0;
            for (int v = max(0, u - h); v <= min(cells - 1, u + h); v++){
                if (v == u || stress::rand_int(0, 2)) continue;
                long double c = stress::rand_int(-1000, 1000) / 250.0L;
                coef[u][v] = c, off += fabsl(c);
            }
            /// Strict diagonal dominance, the library pivots on the diagonal only
            coef[u][u] = it % 4 == 0 && off == 0 ? 1 : (off + stress::rand_int(1, 100) / 10.0L) * (stress::rand_int(0, 1) ? 1 : -1);
            if (stress::rand_int(0, 1)) rhs[u] = stress::rand_int(-1000, 1000) / 8.0L;

            auto [i, j] = cell(u);
            for (auto [v, c] : coef[u]){
                if (v == u && c == 1 && stress::rand_int(0, 1)) continue;  /// the constructor's default diagonal
                auto [k, l] = cell(v);
                gauss->set_matrix(i, j, k, l, c);
            }
            if (rhs[u] != rhs_default || stress::rand_int(0, 1)) gauss->set_rhs(i, j, rhs[u]);
        }

        auto res = gauss->solve();
        assert((int)res.size() == n);
        vector<long double> x(cells);
        for (int i = 0; i < n; i++){
            assert((int)res[i].size() == m);
            for (int j = 0; j < m; j++) x[m * (n - i) - j - 1] = res[i][j];
        }

        long double norm = 1;
        for (auto v : x) norm = max(norm, fabsl(v));
        for (int u = 0; u < cells; u++){
            long double s = 0;
            for (auto [v, c] : coef[u]) s += c * x[v];
            assert(fabsl(s - rhs[u]) <= 1e-12L * norm * (1 + fabsl(coef[u][u])));
        }

        if (cells <= 80){
            vector<vector<__float128>> a(cells, vector<__float128>(cells, 0));
            vector<__float128> b(cells);
            for (int u = 0; u < cells; u++){
                for (auto [v, c] : coef[u]) a[u][v] = c;
                b[u] = rhs[u];
            }
            auto expected = dense_solve(a, b);
            for (int u = 0; u < cells; u++) assert(fabsq(expected[u] - x[u]) <= (__float128)1e-13L * norm);
        }
    }
    operator delete(gauss);
    return 0;
}
