// LINK: -lquadmath
#include "../common.h"
#include <quadmath.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_band_matrix.cpp"
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
    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, it % 10 ? 12 : 400), m = stress::rand_int(1, it % 10 ? 8 : 6);
        if (it % 50 == 0) m = stress::rand_int(31, 60), n = min(n, 30);  /// widest band, kept small: cells times band
        if (stress::rand_int(0, 1)) swap(n, m);  /// so that numbering by column is the narrower band as often as by row

        /// Generated as a band of half width h in its own order, the library must find an order at least as narrow
        bool by_column = stress::rand_int(0, 1);
        int cells = n * m, h = stress::rand_int(0, min(34, cells - 1));
        auto cell = [&](int u){ return by_column ? make_pair(u % n, u / n) : make_pair(u / m, u % m); };

        vector<map<int, long double>> coef(cells);
        vector<long double> rhs(cells, 0);
        vector<bool> untouched(cells, false);
        for (int u = 0; u < cells; u++){
            /// A cell without calls stands for x(u) = 0, other equations may still use it
            if (it % 4 == 0 && stress::rand_int(0, 3) == 0){
                coef[u][u] = 1, untouched[u] = true;
                continue;
            }

            long double off = 0;
            for (int v = max(0, u - h); v <= min(cells - 1, u + h); v++){
                if (v == u || stress::rand_int(0, 2)) continue;
                long double c = stress::rand_int(-1000, 1000) / 250.0L;
                coef[u][v] = c, off += fabsl(c);
            }
            /// Strict diagonal dominance, the library pivots on the diagonal only
            coef[u][u] = (off + stress::rand_int(1, 100) / 10.0L) * (stress::rand_int(0, 1) ? 1 : -1);
            if (stress::rand_int(0, 1)) rhs[u] = stress::rand_int(-1000, 1000) / 8.0L;
        }

        /// sum of a(u, v) x(v) = b(u) is x(u) = b(u) + (1 - a(u, u)) x(u) - sum over v != u of a(u, v) x(v)
        Band gauss(n, m);
        int row_reach = 0, column_reach = 0;
        auto add_term = [&](int i, int j, int k, int l, long double weight){
            gauss.add_term(i, j, k, l, weight);
            row_reach = max(row_reach, abs((i * m + j) - (k * m + l)));
            column_reach = max(column_reach, abs((j * n + i) - (l * n + k)));
        };

        for (int u = 0; u < cells; u++){
            if (untouched[u]) continue;

            auto [i, j] = cell(u);
            for (auto [v, c] : coef[u]){
                long double weight = v == u ? 1 - c : -c;
                if (weight == 0 && stress::rand_int(0, 1)) continue;

                auto [k, l] = cell(v);
                /// Two calls must accumulate
                long double part = stress::rand_int(0, 1) ? stress::rand_int(-1000, 1000) / 64.0L : 0;
                add_term(i, j, k, l, part);
                add_term(i, j, k, l, weight - part);
            }

            long double part = stress::rand_int(0, 1) ? stress::rand_int(-1000, 1000) / 64.0L : 0;
            if (rhs[u] != 0 || stress::rand_int(0, 1)){
                gauss.add_constant(i, j, part);
                gauss.add_constant(i, j, rhs[u] - part);
            }
        }

        auto res = gauss.solve();
        assert(gauss.reach == min(row_reach, column_reach));
        if (row_reach != column_reach) assert(gauss.transposed == (column_reach < row_reach));

        assert((int)res.size() == n);
        vector<long double> x(cells);
        for (int u = 0; u < cells; u++){
            auto [i, j] = cell(u);
            assert((int)res[i].size() == m);
            x[u] = res[i][j];
        }

        long double norm = 1;
        for (auto v : x) norm = max(norm, fabsl(v));
        for (int u = 0; u < cells; u++){
            long double s = 0;
            for (auto [v, c] : coef[u]) s += c * x[v];
            assert(fabsl(s - rhs[u]) <= 1e-12L * norm * (1 + fabsl(coef[u][u])));
            if (untouched[u]) assert(x[u] == 0);
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

    return 0;
}
