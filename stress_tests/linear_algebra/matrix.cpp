#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/matrix.cpp"
#undef main

long long reduce(__int128 x, long long mod){
    x %= mod;
    return (long long)(x < 0 ? x + mod : x);
}

Matrix<> random_matrix(int n, int m, long long mod, int kind){
    Matrix res(n, m, mod);
    for (auto& row : res.a){
        for (auto& x : row){
            if (kind == 0) x = stress::rand_int(0, mod - 1);
            else if (kind == 1) x = stress::rand_int(LLONG_MIN, LLONG_MAX);
            else x = stress::rand_int(-3, 3);
        }
    }

    return res;
}

/// Product against a direct __int128 triple loop
void check_product(const Matrix<>& x, const Matrix<>& y){
    Matrix z = x * y;
    assert(z.n == x.n && z.m == y.m);
    for (int i = 0; i < x.n; i++){
        for (int j = 0; j < y.m; j++){
            __int128 s = 0;
            for (int k = 0; k < x.m; k++) s += (__int128)reduce(x.a[i][k], x.mod) * reduce(y.a[k][j], y.mod);
            assert(z.a[i][j] == reduce(s, x.mod));
        }
    }
}

/// Fibonacci by fast doubling, independent of the matrix code
pair<long long, long long> fib(long long n, long long mod){
    if (n == 0) return {0, 1 % mod};
    auto [a, b] = fib(n >> 1, mod);
    long long c = reduce((__int128)a * reduce(2 * (__int128)b - a, mod), mod);
    long long d = reduce((__int128)a * a + (__int128)b * b, mod);
    if (n & 1) return {d, reduce((__int128)c + d, mod)};
    return {c, d};
}

struct MinPlus{
    static long long zero(){ return LLONG_MAX; }
    static long long one(){ return 0; }
    static long long add(long long x, long long y){ return min(x, y); }
    static long long mul(long long x, long long y){ return x == LLONG_MAX || y == LLONG_MAX ? LLONG_MAX : x + y; }
};

template<typename T>
vector<vector<T>> naive_product(const vector<vector<T>>& x, const vector<vector<T>>& y){
    vector<vector<T>> z(x.size(), vector<T>(y[0].size(), 0));
    for (size_t i = 0; i < x.size(); i++){
        for (size_t j = 0; j < y[0].size(); j++){
            for (size_t k = 0; k < y.size(); k++) z[i][j] += x[i][k] * y[k][j];
        }
    }

    return z;
}

/// Real powers against naive products, with a tolerance relative to the same products over |entries|
void check_real_pow(int n, int e){
    Matrix<long double> x(n, n);
    vector<vector<long double>> slow(n, vector<long double>(n, 0)), bound = slow, magnitude = slow;
    for (int i = 0; i < n; i++){
        slow[i][i] = bound[i][i] = 1;
        for (int j = 0; j < n; j++){
            x[i][j] = stress::rand_int(-1000, 1000) / 128.0L;
            magnitude[i][j] = fabsl(x[i][j]);
        }
    }

    for (int k = 0; k <= e; k++){
        Matrix<long double> p = x.pow(k);
        assert(p.n == n && p.m == n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) assert(fabsl(p[i][j] - slow[i][j]) <= 1e-15L * (1 + bound[i][j]));
        }
        slow = naive_product(slow, x.a);
        bound = naive_product(bound, magnitude);
    }
}

/// Shortest walks with exactly k edges by DP over the edge count, INF where no such walk exists
void check_min_plus(int n, int max_k){
    Matrix<long long, MinPlus> w(n, n);
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++){
            if (stress::rand_int(0, 2)) w[u][v] = stress::rand_int(-1e12, 1e12);
        }
    }

    vector<vector<long long>> dp(n, vector<long long>(n, LLONG_MAX));
    for (int s = 0; s < n; s++) dp[s][s] = 0;
    for (int k = 0; k <= max_k; k++){
        assert(w.pow(k).a == dp);
        vector<vector<long long>> nxt(n, vector<long long>(n, LLONG_MAX));
        for (int s = 0; s < n; s++){
            for (int u = 0; u < n; u++){
                if (dp[s][u] == LLONG_MAX) continue;
                for (int v = 0; v < n; v++){
                    if (w[u][v] != LLONG_MAX) nxt[s][v] = min(nxt[s][v], dp[s][u] + w[u][v]);
                }
            }
        }
        dp = nxt;
    }
}

int main(){
    const long long MAXMOD = (1LL << 31) - 1;

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long mod = it % 4 == 0 ? MAXMOD : (it % 4 == 1 ? stress::rand_int(1, 10) : stress::rand_int(1, MAXMOD));
        int n = stress::rand_int(1, 7), m = stress::rand_int(1, 7), p = stress::rand_int(1, 7);
        check_product(random_matrix(n, m, mod, it % 3), random_matrix(m, p, mod, (it + 1) % 3));
    }

    /// multiply must overwrite whatever an earlier product left in res
    for (long long it = 0; it < stress::scaled(1000); it++){
        long long mod = it % 2 ? MAXMOD : stress::rand_int(1, 20);
        int n = stress::rand_int(1, 6), m = stress::rand_int(1, 6), p = stress::rand_int(1, 6);
        Matrix x = random_matrix(n, m, mod, it % 3), y = random_matrix(m, p, mod, 0);
        Matrix res = random_matrix(n, p, mod, 1);
        Matrix<>::multiply(x, y, res);
        assert(res.a == (x * y).a);

        bool reduced = true;
        for (auto& row : x.a){
            for (long long v : row) reduced &= 0 <= v && v < mod;
        }
        assert(x.is_normalized() == reduced);
    }

    for (long long it = 0; it < stress::scaled(30); it++){
        long long mod = it % 2 ? MAXMOD : stress::rand_int(1, MAXMOD);
        int n = stress::rand_int(30, 70);
        check_product(random_matrix(n, n, mod, it % 3), random_matrix(n, n, mod, 0));
    }

    for (long long it = 0; it < stress::scaled(500); it++){
        long long mod = it % 3 ? MAXMOD : (it % 9 ? stress::rand_int(1, 50) : 1);
        int n = stress::rand_int(1, 5);
        Matrix x = random_matrix(n, n, mod, it % 3);
        Matrix slow(n, n, mod);  /// built by hand so a wrong Matrix::identity cannot hide in both sides
        for (int i = 0; i < n; i++) slow.a[i][i] = reduce(1, mod);
        for (int e = 0; e <= 40; e++){
            assert(x.pow(e).a == slow.a);
            slow = slow * x;
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long mod = it % 3 ? MAXMOD : stress::rand_int(1, MAXMOD);
        long long n = it % 5 == 0 ? LLONG_MAX - stress::rand_int(0, 1000) : stress::rand_int(0, LLONG_MAX);
        Matrix q(2, 2, mod);
        q.a = {{1, 1}, {1, 0}};
        assert(q.pow(n)[0][1] == fib(n, mod).first);
    }

    for (long long it = 0; it < stress::scaled(500); it++){
        int n = stress::rand_int(1, 8), m = stress::rand_int(1, 8);
        Matrix x = random_matrix(n, m, MAXMOD, 0);
        Matrix t = x.transpose(), r = x.rotate();
        assert(t.n == m && t.m == n && r.n == m && r.m == n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                assert(t.a[j][i] == x.a[i][j]);
                assert(r.a[j][n - 1 - i] == x.a[i][j]);
            }
        }
        assert(x.rotate().rotate().rotate().rotate().a == x.a);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 6), m = stress::rand_int(1, 6), p = stress::rand_int(1, 6);
        Matrix<long double> x(n, m), y(m, p);
        for (auto& row : x.a){
            for (auto& v : row) v = stress::rand_int(-1e6, 1e6) / 1024.0L;
        }
        for (auto& row : y.a){
            for (auto& v : row) v = stress::rand_int(-1e6, 1e6) / 1024.0L;
        }
        /// multiples of 2^-10 below 2^10, so every product and sum is exact
        assert((x * y).a == naive_product(x.a, y.a));
        check_real_pow(stress::rand_int(1, 6), 20);
    }

    for (long long it = 0; it < stress::scaled(300); it++) check_min_plus(stress::rand_int(1, 6), 12);

    /// No modulus: long long stays exact and unsigned long long wraps modulo 2^64, never reduced by a stale mod
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 6), m = stress::rand_int(1, 6), p = stress::rand_int(1, 6);
        Matrix<long long> x(n, m), y(m, p);
        Matrix<unsigned long long> ux(n, m), uy(m, p);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) x[i][j] = stress::rand_int(-1e9, 1e9), ux[i][j] = stress::rng()();
        }
        for (int i = 0; i < m; i++){
            for (int j = 0; j < p; j++) y[i][j] = stress::rand_int(-1e9, 1e9), uy[i][j] = stress::rng()();
        }
        assert(x.is_normalized());
        assert(x.normalized().a == x.a);
        assert((x * y).a == naive_product(x.a, y.a));
        assert((ux * uy).a == naive_product(ux.a, uy.a));
        assert(Matrix<long long>::identity(n).pow(5).a == Matrix<long long>::identity(n).a);
    }

    stress::expect_assert("linear_algebra/matrix.cpp", []{ Matrix<long double>(1, 1, 7); });
    stress::expect_assert("linear_algebra/matrix.cpp", []{ Matrix<long long, MinPlus>(1, 1, 7); });
    stress::expect_assert("linear_algebra/matrix.cpp", []{ Matrix<>(1, 1, 1LL << 31); });

    return 0;
}
