/***
 *
 * Matrix
 * Modular matrix multiplication and exponentiation, plus transpose and rotation
 *
 * Complexity: O(n * m * p) per multiplication, O(n^3 log e) for pow(e)
 *
 * Matrix a(n, m, mod): n x m zero matrix modulo mod, 1 <= mod < 2^31, Matrix::identity(n, mod)
 * a[i][j] = value: any long long, negative values are reduced modulo mod when multiplying
 * a * b: needs a.m == b.n and the same mod
 * a.pow(e): a^e for a square matrix and e >= 0, a.pow(0) is the identity
 * a.transpose(), a.rotate(): rotate turns the matrix 90 degrees clockwise
 *
 * Linear recurrences: f(n) = 2 f(n - 1) + 3 f(n - 2) uses M = {{2, 3}, {1, 0}},
 * then (f(n + 1), f(n)) = M^n * (f(1), f(0))
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Matrix{
    int n, m;
    long long mod;
    vector<vector<long long>> a;

    Matrix(int n, int m, long long mod) : n(n), m(m), mod(mod), a(n, vector<long long>(m, 0)){
        assert(1 <= mod && mod < (1LL << 31));
    }

    static Matrix identity(int n, long long mod){
        Matrix res(n, n, mod);
        for (int i = 0; i < n; i++) res.a[i][i] = 1 % mod;
        return res;
    }

    vector<long long>& operator[](int i){
        return a[i];
    }

    const vector<long long>& operator[](int i) const{
        return a[i];
    }

    Matrix normalized() const{
        Matrix res = *this;
        for (auto& row : res.a){
            for (auto& x : row){
                x %= mod;
                if (x < 0) x += mod;
            }
        }

        return res;
    }

    Matrix operator*(const Matrix& other) const{
        assert(m == other.n && mod == other.mod);

        Matrix x = normalized(), y = other.normalized(), res(n, other.m, mod);
        const long long limit = mod * mod;  /// sums stay below 2 * mod^2 < 2^63
        vector<long long> acc(other.m);

        for (int i = 0; i < n; i++){
            fill(acc.begin(), acc.end(), 0);
            for (int k = 0; k < m; k++){
                long long v = x.a[i][k];
                if (!v) continue;
                for (int j = 0; j < other.m; j++){
                    acc[j] += v * y.a[k][j];
                    if (acc[j] >= limit) acc[j] -= limit;
                }
            }
            for (int j = 0; j < other.m; j++) res.a[i][j] = acc[j] % mod;
        }

        return res;
    }

    Matrix pow(long long e) const{
        assert(n == m && e >= 0);

        Matrix res = identity(n, mod), base = normalized();
        for (; e; e >>= 1){
            if (e & 1) res = res * base;
            base = base * base;
        }

        return res;
    }

    Matrix transpose() const{
        Matrix res(m, n, mod);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) res.a[j][i] = a[i][j];
        }

        return res;
    }

    Matrix rotate() const{
        Matrix res = transpose();
        for (auto& row : res.a) reverse(row.begin(), row.end());
        return res;
    }
};

const long long MOD = 10000007;

int R;

void fill_column(Matrix& t, int in, int row, int filled, int out){
    if (row == R){
        t[in][out]++;
        return;
    }
    if (filled >> row & 1) return fill_column(t, in, row + 1, filled, out);

    fill_column(t, in, row + 1, filled, out | 1 << row);
    if (row + 1 < R && !(filled >> (row + 1) & 1)) fill_column(t, in, row + 2, filled, out);
}

Matrix transfer(int blocked){
    Matrix t(1 << R, 1 << R, MOD);
    for (int in = 0; in < (1 << R); in++){
        if (!(in & blocked)) fill_column(t, in, 0, in | blocked, 0);
    }

    return t;
}

int main(){
    int C, n, tc = 0;

    while (scanf("%d %d %d", &R, &C, &n) == 3 && R){
        map<int, int> blocked;
        for (int i = 0; i < n; i++){
            int r, c;
            scanf("%d %d", &r, &c);
            blocked[c] |= 1 << r;
        }

        vector<Matrix> free_pow = {transfer(0)};
        while ((1 << free_pow.size()) <= C) free_pow.push_back(free_pow.back() * free_pow.back());

        Matrix v(1, 1 << R, MOD);
        v[0][0] = 1;
        auto advance = [&](int columns){
            for (int k = 0; columns; k++, columns >>= 1){
                if (columns & 1) v = v * free_pow[k];
            }
        };

        int done = 0;
        for (auto e : blocked){
            int c = e.first, mask = e.second;
            advance(c - done);
            v = v * transfer(mask);
            done = c + 1;
        }
        advance(C - done);
        printf("Case %d: %lld\n", ++tc, v[0][0]);
    }

    return 0;
}
