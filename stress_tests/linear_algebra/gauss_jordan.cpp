#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/gauss_jordan.cpp"
#undef main

/// Exact rank via elimination modulo a large prime, it matches the rational rank unless the prime divides a minor
int rank_mod(const vector<vector<long long>>& a, int cols, long long p){
    int rows = a.size(), rank = 0;
    vector<vector<long long>> b(rows, vector<long long>(cols));
    for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) b[i][j] = ((a[i][j] % p) + p) % p;
    auto mul = [&](long long x, long long y){ return (long long)((__int128)x * y % p); };
    auto inv = [&](long long x){ long long r = 1; for (long long e = p - 2; e; e >>= 1, x = mul(x, x)) if (e & 1) r = mul(r, x); return r; };

    for (int c = 0; c < cols && rank < rows; c++){
        int piv = rank;
        while (piv < rows && b[piv][c] == 0) piv++;
        if (piv == rows) continue;
        swap(b[piv], b[rank]);
        long long f = inv(b[rank][c]);
        for (int r = rank + 1; r < rows; r++){
            long long t = mul(b[r][c], f);
            for (int k = c; k < cols; k++) b[r][k] = (b[r][k] - mul(t, b[rank][k]) + p) % p;
        }
        rank++;
    }
    return rank;
}

int exact_rank(const vector<vector<long long>>& a, int cols){
    int r1 = rank_mod(a, cols, (1LL << 61) - 1), r2 = rank_mod(a, cols, 1000000000000000003LL);
    assert(r1 == r2);
    return r1;
}

int main(){
    /// Consistent (row 5 = row 1 + row 4) but rejected by the old 1e-12 tolerance, found by this test at 20x scale
    vector<vector<double>> regression = {{461376, -874886, 67590, 302992, 246673}, {257837, 571397, 503158, -220469, 555968},
                                         {-906930, 685565, -202946, 547555, 179852}, {-16614, 203028, -207928, -632214, 573819},
                                         {444762, -671858, -140338, -329222, 820492}};
    vector<double> solution;
    assert(gauss(regression, solution) == 0);

    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(1, 5), m = stress::rand_int(1, 5);
        long long big = stress::rand_int(0, 2) ? 10 : 1000000;
        vector<vector<long long>> a(n, vector<long long>(m + 1));
        for (int i = 0; i < n; i++){
            int kind = i ? stress::rand_int(0, 3) : 0;
            for (int j = 0; j <= m; j++){
                if (kind == 0 || kind == 1) a[i][j] = stress::rand_int(-big, big);         /// independent row
                else if (kind == 2) a[i][j] = a[0][j] * 3;                                  /// dependent row
                else a[i][j] = a[0][j] + a[i - 1][j];                                       /// combination of rows
            }
            if (kind >= 2 && stress::rand_int(0, 2) == 0) a[i][m] += stress::rand_int(1, 5);  /// breaks consistency
        }

        vector<vector<double>> eq(n, vector<double>(m + 1));
        for (int i = 0; i < n; i++) for (int j = 0; j <= m; j++) eq[i][j] = a[i][j];
        int rank = exact_rank(a, m), expected = rank == exact_rank(a, m + 1) ? m - rank : -1;

        vector<double> res;
        int got = gauss(eq, res);
        assert(got == expected);
        if (got >= 0){
            for (int i = 0; i < n; i++){
                double lhs = 0, scale = fabs(eq[i][m]);
                for (int j = 0; j < m; j++) lhs += eq[i][j] * res[j], scale = max(scale, fabs(eq[i][j] * res[j]));
                assert(fabs(lhs - eq[i][m]) <= 1e-9 * max(1.0, scale));
            }
        }
    }
    return 0;
}
