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
        int n = stress::rand_int(1, it % 10 ? 12 : 3000);
        vector<equation> ar(n);
        for (int i = 0; i < n; i++){
            /// Coefficients outside the matrix must be ignored, so they get junk
            long double l = stress::rand_int(-1000, 1000) / 100.0L, r = stress::rand_int(-1000, 1000) / 100.0L;
            long double p = (fabsl(i ? l : 0) + fabsl(i + 1 < n ? r : 0) + stress::rand_int(1, 1000) / 100.0L) * (stress::rand_int(0, 1) ? 1 : -1);
            ar[i] = equation(l, p, r, stress::rand_int(-100000, 100000) / 64.0L);
        }

        auto x = thomas_algorithm(n, ar);
        assert((int)x.size() == n);

        long double norm = 1;
        for (auto v : x) norm = max(norm, fabsl(v));
        for (int i = 0; i < n; i++){
            long double s = ar[i].p * x[i] + (i ? ar[i].l * x[i - 1] : 0) + (i + 1 < n ? ar[i].r * x[i + 1] : 0);
            assert(fabsl(s - ar[i].rhs) <= 1e-17L * norm * (fabsl(ar[i].p) + fabsl(ar[i].l) + fabsl(ar[i].r)));
        }

        if (n <= 60){
            vector<vector<__float128>> a(n, vector<__float128>(n, 0));
            vector<__float128> b(n);
            for (int i = 0; i < n; i++){
                a[i][i] = ar[i].p;
                if (i) a[i][i - 1] = ar[i].l;
                if (i + 1 < n) a[i][i + 1] = ar[i].r;
                b[i] = ar[i].rhs;
            }
            auto expected = dense_solve(a, b);
            for (int i = 0; i < n; i++) assert(fabsq(expected[i] - x[i]) <= (__float128)1e-17L * norm);
        }
    }

    return 0;
}
