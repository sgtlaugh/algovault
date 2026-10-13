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

int main(){
    const long long MOD = 10007;
    vector<int> cols;
    for (int mask = 0; mask < 64; mask++){
        if (__builtin_popcount(mask) == 2) cols.push_back(mask);
    }

    /// Knights one column apart attack across 2 rows, two columns apart across 1 row
    auto apart = [](int a, int b, int d){ return !((a << d) & b) && !((b << d) & a); };

    /// A state is the pair of the last two columns, the next column must suit both
    vector<pair<int, int>> states;
    for (int a : cols){
        for (int b : cols){
            if (apart(a, b, 2)) states.push_back({a, b});
        }
    }

    int s = states.size();
    Matrix step(s, s, MOD), start(1, s, MOD);
    for (int i = 0; i < s; i++){
        start[0][i] = 1;
        for (int j = 0; j < s; j++){
            auto [a, b] = states[i];
            auto [c, d] = states[j];
            if (b == c && apart(a, d, 1)) step[i][j] = 1;
        }
    }

    /// Squarings shared by all test cases, each query only multiplies the row vector through
    vector<Matrix> powers = {step};
    for (int k = 1; k < 30; k++) powers.push_back(powers.back() * powers.back());

    int t;
    scanf("%d", &t);
    while (t--){
        int n;
        scanf("%d", &n);
        if (n == 1){
            printf("%d\n", (int)cols.size());
            continue;
        }

        Matrix v = start;
        for (int k = 0; k < 30; k++){
            if ((n - 2) >> k & 1) v = v * powers[k];
        }
        long long res = 0;
        for (int i = 0; i < s; i++) res += v[0][i];
        printf("%lld\n", res % MOD);
    }
    return 0;
}
