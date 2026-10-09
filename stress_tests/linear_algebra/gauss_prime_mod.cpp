#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/gauss_prime_mod.cpp"
#undef main

bool satisfies(const vector<vector<int>>& equations, const vector<int>& x, int mod){
    int m = x.size();
    for (const auto& e : equations){
        long long v = 0;
        for (int j = 0; j < m; j++) v = (v + (long long)e[j] * x[j]) % mod;
        if (v != e[m]) return false;
    }
    return true;
}

vector<vector<int>> homogeneous(vector<vector<int>> equations){
    for (auto& e : equations) e.back() = 0;
    return equations;
}

/// Rank by forward elimination with first-nonzero pivots, written independently of the library
int rank_ref(vector<vector<long long>> a, int cols, long long p){
    int rows = a.size(), rank = 0;
    auto power = [&](long long x, long long e){ long long r = 1; for (; e; e >>= 1, x = x * x % p) if (e & 1) r = r * x % p; return r; };

    for (int c = 0; c < cols && rank < rows; c++){
        int piv = rank;
        while (piv < rows && a[piv][c] % p == 0) piv++;
        if (piv == rows) continue;
        swap(a[piv], a[rank]);
        long long f = power((a[rank][c] % p + p) % p, p - 2);
        for (int r = rank + 1; r < rows; r++){
            long long t = (a[r][c] % p + p) % p * f % p;
            for (int k = c; k < cols; k++) a[r][k] = ((a[r][k] - t * (a[rank][k] % p)) % p + p) % p;
        }
        rank++;
    }

    return rank;
}

vector<vector<long long>> widen(const vector<vector<int>>& a){
    vector<vector<long long>> res;
    for (const auto& row : a) res.emplace_back(row.begin(), row.end());
    return res;
}

bool is_identity_product(const vector<vector<int>>& x, const vector<vector<int>>& y, int mod){
    int n = x.size();
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            long long v = 0;
            for (int k = 0; k < n; k++) v = (v + (long long)((x[i][k] % mod + mod) % mod) * ((y[k][j] % mod + mod) % mod)) % mod;
            if (v != (i == j)) return false;
        }
    }
    return true;
}

/// res plus every combination of the basis must be exactly the solution set, each solution reached once
void check_solution_space(const vector<vector<int>>& equations, const vector<int>& res, const vector<vector<int>>& basis, int mod, long long solutions){
    int m = res.size(), f = basis.size();
    long long combos = 1;
    for (int k = 0; k < f; k++) combos *= mod;
    assert(combos == solutions);

    set<vector<int>> seen;
    vector<int> c(f, 0);
    for (long long code = 0; code < combos; code++){
        long long t = code;
        for (int k = 0; k < f; k++) c[k] = t % mod, t /= mod;
        vector<int> x = res;
        for (int k = 0; k < f; k++){
            for (int j = 0; j < m; j++) x[j] = (x[j] + (long long)c[k] * basis[k][j]) % mod;
        }
        assert(satisfies(equations, x, mod));
        seen.insert(x);
    }
    assert((long long)seen.size() == combos);
}

int main(){
    vector<int> res;
    assert(gauss({{1, 8}}, res, 7) == 0 && res[0] == 1);  /// unreduced RHS
    assert(gauss({{7, 0}}, res, 7) == 1 && res[0] == 0);  /// coefficient 7 is zero mod 7
    assert(gauss({{7, 3}}, res, 7) == -1);
    assert(gauss({{-1, -3}, {-2, 1}}, res, 5) == -1);  /// -x = -3 and -2x = 1 imply x = 3 and x = 2
    assert(gauss({{-1, 0, -3}, {0, -2, 1}}, res, 5) == 0 && res[0] == 3 && res[1] == 2);
    assert(gauss({{INT_MIN, INT_MAX, 0}}, res, 1000000007) == 1 && satisfies({{(INT_MIN % 1000000007) + 1000000007, INT_MAX % 1000000007, 0}}, res, 1000000007));

    const int max_prime = 2147483647;
    vector<vector<int>> inv;
    assert(matrix_inverse({{max_prime - 1, 0}, {0, 2}}, inv, max_prime) && (inv == vector<vector<int>>{{max_prime - 1, 0}, {0, max_prime / 2 + 1}}));
    assert(matrix_rank({{INT_MAX, 1}, {0, 1}}, max_prime) == 1);

    const int small_primes[] = {2, 3, 5, 7};
    for (long long it = 0; it < stress::scaled(3000); it++){
        int mod = small_primes[stress::rand_int(0, 3)];
        int m = stress::rand_int(1, mod == 2 ? 10 : mod == 3 ? 6 : 4), n = stress::rand_int(1, 7);
        vector<vector<int>> equations(n, vector<int>(m + 1));
        int density = stress::rand_int(0, 2);
        for (auto& e : equations) for (auto& x : e) x = stress::rand_int(0, density) ? stress::rand_int(0, mod - 1) : 0;
        if (n >= 2 && stress::rand_int(0, 1)){
            int r = stress::rand_int(1, n - 1), c = stress::rand_int(1, mod - 1);
            for (int j = 0; j <= m; j++) equations[r][j] = (equations[r][j] + c * equations[0][j]) % mod;  /// a dependent row
        }

        long long solutions = 0, kernel = 0, total = 1;
        for (int j = 0; j < m; j++) total *= mod;
        vector<vector<int>> zero_rhs = homogeneous(equations);
        vector<int> x(m);
        for (long long code = 0; code < total; code++){
            long long c = code;
            for (int j = 0; j < m; j++) x[j] = c % mod, c /= mod;
            solutions += satisfies(equations, x, mod);
            kernel += satisfies(zero_rhs, x, mod);
        }

        vector<vector<int>> unreduced = equations;  /// any int is congruent to its reduction
        if (stress::rand_int(0, 1)){
            for (auto& e : unreduced) for (auto& v : e) v += mod * stress::rand_int(-3, 3);
        }

        vector<int> res = {42, 42};  /// stale contents must be replaced
        vector<vector<int>> basis = {{1}};
        int f_var = gauss(unreduced, res, basis, mod);
        check_solution_space(zero_rhs, vector<int>(m, 0), basis, mod, kernel);
        if (!solutions) assert(f_var == -1);
        else{
            long long expected = 1;
            for (int j = 0; j < f_var; j++) expected *= mod;
            assert(f_var >= 0 && solutions == expected);
            assert((int)res.size() == m && satisfies(equations, res, mod));
            for (int v : res) assert(0 <= v && v < mod);
            for (const auto& b : basis) for (int v : b) assert(0 <= v && v < mod);
            check_solution_space(equations, res, basis, mod, solutions);
            assert(gauss(unreduced, res, mod) == f_var && satisfies(equations, res, mod));
        }

        /// The kernel has mod^(m - rank) vectors
        vector<vector<int>> coefficients = unreduced;
        for (auto& e : coefficients) e.pop_back();
        long long size = 1;
        int rank = matrix_rank(coefficients, mod);
        for (int j = 0; j < m - rank; j++) size *= mod;
        assert(size == kernel);
    }

    /// Small square matrices: invertible exactly when the kernel is trivial, checked by enumeration
    for (long long it = 0; it < stress::scaled(3000); it++){
        int mod = small_primes[stress::rand_int(0, 3)];
        int n = stress::rand_int(0, mod == 2 ? 6 : mod == 3 ? 4 : 3), density = stress::rand_int(0, 3);
        vector<vector<int>> a(n, vector<int>(n));
        for (auto& row : a) for (auto& v : row) v = stress::rand_int(0, density) ? stress::rand_int(-20, 20) : 0;

        long long total = 1, kernel = 0;
        for (int j = 0; j < n; j++) total *= mod;
        for (long long code = 0; code < total; code++){
            long long c = code;
            bool zero = true;
            for (int i = 0; i < n && zero; i++){
                long long v = 0;
                c = code;
                for (int j = 0; j < n; j++, c /= mod) v += (long long)a[i][j] * (c % mod);
                zero = v % mod == 0;
            }
            kernel += zero;
        }

        vector<vector<int>> inv = {{3}};
        bool ok = matrix_inverse(a, inv, mod);
        assert(ok == (kernel == 1));
        if (ok){
            assert((int)inv.size() == n && is_identity_product(a, inv, mod) && is_identity_product(inv, a, mod));
            for (const auto& row : inv) for (int v : row) assert(0 <= v && v < mod);
        }
    }

    /// A large prime: rank r systems built as products of n x r and r x m random factors, full rank r almost surely
    const int mod = 1000000007;
    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(1, it % 10 ? 10 : 60), m = stress::rand_int(1, it % 10 ? 10 : 60), r = stress::rand_int(0, min(n, m));
        vector<vector<long long>> left(n, vector<long long>(r)), right(r, vector<long long>(m));
        for (auto& row : left) for (auto& v : row) v = stress::rand_int(0, mod - 1);
        for (auto& row : right) for (auto& v : row) v = stress::rand_int(0, mod - 1);
        vector<int> planted(m);
        for (auto& v : planted) v = stress::rand_int(0, mod - 1);

        vector<vector<int>> equations(n, vector<int>(m + 1));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                long long v = 0;
                for (int k = 0; k < r; k++) v = (v + left[i][k] * right[k][j]) % mod;
                equations[i][j] = v;
            }
            long long v = 0;
            for (int j = 0; j < m; j++) v = (v + (long long)equations[i][j] * planted[j]) % mod;
            equations[i][m] = v;
        }

        if (r < n){
            /// The other n - 1 >= r rows span the row space almost surely, so shifting any RHS makes the system inconsistent
            vector<vector<int>> probe = equations;
            int row = stress::rand_int(0, n - 1);
            probe[row][m] = (probe[row][m] + stress::rand_int(1, mod - 1)) % mod;
            vector<int> res;
            assert(gauss(probe, res, mod) == -1);
        }

        vector<int> res;
        int f_var = gauss(equations, res, mod);
        assert(f_var == m - r && satisfies(equations, res, mod));

        vector<vector<int>> basis;
        assert(gauss(equations, res, basis, mod) == m - r && (int)basis.size() == m - r);
        assert(rank_ref(widen(basis), m, mod) == m - r);
        for (int t = 0; t < 5; t++){
            vector<int> x = res;
            for (const auto& b : basis){
                long long c = stress::rand_int(0, mod - 1);
                for (int j = 0; j < m; j++) x[j] = (x[j] + c * b[j]) % mod;
            }
            assert(satisfies(equations, x, mod));
        }

        vector<vector<int>> coefficients = equations;
        for (auto& e : coefficients) e.pop_back();
        int rank = matrix_rank(coefficients, mod);
        assert(rank == rank_ref(widen(coefficients), m, mod) && rank == r);

        /// Square matrices: random ones are invertible almost surely, the n x r times r x n products with r < n never are
        vector<vector<int>> square(n, vector<int>(n));
        for (auto& row : square) for (auto& v : row) v = stress::rand_int(INT_MIN, INT_MAX);
        if (r < n){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    long long v = 0;
                    for (int k = 0; k < r; k++) v = (v + left[i][k] * left[j][k]) % mod;
                    square[i][j] = v;
                }
            }
        }

        vector<vector<int>> inv;
        bool ok = matrix_inverse(square, inv, mod);
        assert(ok == (rank_ref(widen(square), n, mod) == n));
        if (ok) assert(is_identity_product(square, inv, mod) && is_identity_product(inv, square, mod));
    }

    return 0;
}
