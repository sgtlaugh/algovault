#include "common.h"

#define main library_main
#include "../code_library/matrix.cpp"
#undef main

long long reduce(__int128 x, long long mod){
    x %= mod;
    return (long long)(x < 0 ? x + mod : x);
}

Matrix random_matrix(int n, int m, long long mod, int kind){
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
void check_product(const Matrix& x, const Matrix& y){
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

int main(){
    const long long MAXMOD = (1LL << 31) - 1;

    for (long long it = 0; it < stress::scaled(3000); it++){
        long long mod = it % 4 == 0 ? MAXMOD : (it % 4 == 1 ? stress::rand_int(1, 10) : stress::rand_int(1, MAXMOD));
        int n = stress::rand_int(1, 7), m = stress::rand_int(1, 7), p = stress::rand_int(1, 7);
        check_product(random_matrix(n, m, mod, it % 3), random_matrix(m, p, mod, (it + 1) % 3));
    }
    for (long long it = 0; it < stress::scaled(30); it++){
        long long mod = it % 2 ? MAXMOD : stress::rand_int(1, MAXMOD);
        int n = stress::rand_int(30, 70);
        check_product(random_matrix(n, n, mod, it % 3), random_matrix(n, n, mod, 0));
    }

    for (long long it = 0; it < stress::scaled(500); it++){
        long long mod = it % 3 ? MAXMOD : stress::rand_int(1, 50);
        int n = stress::rand_int(1, 5);
        Matrix x = random_matrix(n, n, mod, it % 3);
        Matrix slow = Matrix::identity(n, mod);
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

    return 0;
}
