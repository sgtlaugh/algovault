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

const long long MOD = 34830;

int len;
long long k, dp[12][3][3];

/// Counts the ways to color chairs i..len - 1 when u and v are the colors of the two chairs before: 0 or 1 for the color
/// of chair 0 or chair 1, 2 for any other, since only those two matter when the circle closes
long long ways(int i, int u, int v){
    if (i == len) return (u == 0 || v == 0 || v == 1) ? 0 : 1;
    if (dp[i][u][v] != -1) return dp[i][u][v];

    long long res = 0;
    if ((u == 0 && v == 1) || (u == 1 && v == 0)) res += (k - 2) * ways(i + 1, v, 2);
    if ((u == 0 && v == 2) || (u == 1 && v == 2)) res += ways(i + 1, v, u ^ 1) + (k - 3) * ways(i + 1, v, 2);
    if (u == 2 && v != 2) res += ways(i + 1, v, v ^ 1) + (k - 3) * ways(i + 1, v, 2);
    if (u == 2 && v == 2) res += ways(i + 1, v, 0) + ways(i + 1, v, 1) + (k - 4) * ways(i + 1, v, 2);
    return dp[i][u][v] = res % MOD;
}

long long special_case(long long n){
    if (n == 1) return k;
    if (k == 1) return 0;
    if (n == 2 && k == 2) return 2;
    if (k == 2) return 0;
    if (n == 2 && k == 3) return 6;
    if (k == 3) return n % 3 == 0 ? 6 : 0;
    return -1;
}

long long small_case(int n){
    if (special_case(n) != -1) return special_case(n);

    /// Up to 5 chairs every pair lies within some three consecutive ones, so all colors differ
    if (n <= 5){
        long long res = 1;
        for (int j = 0; j < n; j++) res = res * ((k - j) % MOD) % MOD;
        return res;
    }

    len = n;
    memset(dp, -1, sizeof(dp));
    long long x = 1, y = 1;
    for (int j = 0; j < 3; j++) y = y * ((k - j) % MOD) % MOD;
    for (int j = 0; j < 4; j++) x = x * ((k - j) % MOD) % MOD;
    return (x * ways(4, 2, 2) + y * ways(4, 2, 0)) % MOD;
}

int main(){
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        long long n;
        scanf("%lld %lld", &n, &k);

        long long res = special_case(n);
        if (res == -1 && n <= 7) res = small_case(n);
        else if (res == -1){
            /// f(n) = (k - 3) f(n - 1) + f(n - 2) + ((k - 2)^2 - (k - 3)) f(n - 3) - (k - 2)^2 f(n - 5) for n > 7
            long long c = (k - 2) % MOD * ((k - 2) % MOD) % MOD;
            Matrix step(5, 5, MOD);
            step[0][0] = k - 3, step[0][1] = 1, step[0][2] = c - (k - 3), step[0][4] = -c;
            for (int i = 1; i < 5; i++) step[i][i - 1] = 1;

            Matrix power = step.pow(n - 7);
            res = 0;
            for (int i = 0; i < 5; i++) res = (res + power[0][i] * small_case(7 - i)) % MOD;
        }
        printf("Case %d: %lld\n", tc, res % MOD);
    }
    return 0;
}
