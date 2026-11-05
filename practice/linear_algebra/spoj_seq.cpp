/***
 *
 * Matrix
 * Matrix products and powers modulo mod or over any semiring, plus transpose and rotation
 *
 * Complexity:
 *   - multiplication: O(n m p) for an n x m times an m x p matrix
 *   - pow(e): O(n^3 log e)
 *
 * Matrix a(n, m, mod): n x m zero matrix of long long modulo mod, 1 <= mod < 2^31
 * Matrix<T, Ops = PlusTimes<T>> a(n, m): n x m matrix of Ops::zero() with no modulus, like Matrix<long double> a(n, m)
 * Matrix<T, Ops>::identity(n, mod = 0): the n x n identity, written Matrix<>::identity(n, mod) when modular
 * a[i][j] = value: any T; with a modulus any long long, reduced modulo mod when multiplying
 * a * b: the product, needs a.m == b.n and the same mod
 * a.is_normalized(): whether every entry is already in [0, mod), always true without a modulus
 * a.normalized(): a copy with every entry in [0, mod), an unchanged copy without a modulus
 * a.pow(e): a^e for a square matrix and e >= 0, a.pow(0) is the identity
 * a.rotate(): the matrix turned 90 degrees clockwise
 * a.transpose(): the m x n transpose
 *
 * Notes:
 *   - Recurrence f(n) = 2 f(n - 1) + 3 f(n - 2): with M = {{2, 3}, {1, 0}}, (f(n + 1), f(n)) = M^n * (f(1), f(0))
 *   - Ops holds static zero(), one(), add(x, y), mul(x, y) of a semiring, so zero must absorb in mul
 *   - Min-plus: zero = INF, one = 0, add = min, mul = + capped at INF; a.pow(k)[u][v] is the shortest k-edge walk
 *   - Only the default Matrix<long long> takes a modulus, any other T or Ops asserts mod == 0
 *   - Without a modulus nothing is reduced: long long can overflow, unsigned long long wraps modulo 2^64
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T>
struct PlusTimes{
    static T zero(){ return 0; }
    static T one(){ return 1; }
    static T add(T x, T y){ return x + y; }
    static T mul(T x, T y){ return x * y; }
};

template<typename T = long long, typename Ops = PlusTimes<T>>
struct Matrix{
    static constexpr bool MODULAR = is_same<T, long long>::value && is_same<Ops, PlusTimes<T>>::value;

    int n, m;
    long long mod;
    vector<vector<T>> a;

    Matrix(int n, int m, long long mod = 0) : n(n), m(m), mod(mod), a(n, vector<T>(m, Ops::zero())){
        assert(mod == 0 || (MODULAR && 1 <= mod && mod < (1LL << 31)));
    }

    static Matrix identity(int n, long long mod = 0){
        Matrix res(n, n, mod);
        for (int i = 0; i < n; i++) res.a[i][i] = Ops::one();
        return res.normalized();
    }

    vector<T>& operator[](int i){
        return a[i];
    }

    const vector<T>& operator[](int i) const{
        return a[i];
    }

    Matrix operator*(const Matrix& other) const{
        assert(m == other.n && mod == other.mod);

        Matrix res(n, other.m, mod);
        if (other.is_normalized()) multiply(*this, other, res);
        else multiply(*this, other.normalized(), res);
        return res;
    }

    bool is_normalized() const{
        if constexpr (MODULAR){
            if (!mod) return true;
            for (const auto& row : a){
                for (long long x : row){
                    if (x < 0 || x >= mod) return false;
                }
            }
        }

        return true;
    }

    Matrix normalized() const{
        Matrix res = *this;

        if constexpr (MODULAR){
            if (!mod) return res;
            for (auto& row : res.a){
                for (auto& x : row){
                    x %= mod;
                    if (x < 0) x += mod;
                }
            }
        }

        return res;
    }

    Matrix pow(long long e) const{
        assert(n == m && e >= 0);

        /// tmp is the one scratch buffer, so the loop allocates nothing
        Matrix res = identity(n, mod), base = normalized(), tmp(n, n, mod);
        for (; e; e >>= 1){
            if (e & 1){
                multiply(res, base, tmp);
                swap(res, tmp);
            }
            multiply(base, base, tmp);
            swap(base, tmp);
        }

        return res;
    }

    Matrix rotate() const{
        Matrix res = transpose();
        for (auto& row : res.a) reverse(row.begin(), row.end());
        return res;
    }

    Matrix transpose() const{
        Matrix res(m, n, mod);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++) res.a[j][i] = a[i][j];
        }

        return res;
    }

    /// res = x * y in res's existing rows; y must be normalized and res must not alias x or y
    static void multiply(const Matrix& x, const Matrix& y, Matrix& res){
        if constexpr (MODULAR){
            if (x.mod) return multiply_mod(x, y, res);
        }

        for (int i = 0; i < x.n; i++){
            vector<T>& acc = res.a[i];
            fill(acc.begin(), acc.end(), Ops::zero());
            for (int k = 0; k < x.m; k++){
                const T v = x.a[i][k];
                const vector<T>& row = y.a[k];
                for (int j = 0; j < y.m; j++) acc[j] = Ops::add(acc[j], Ops::mul(v, row[j]));
            }
        }
    }

    static void multiply_mod(const Matrix& x, const Matrix& y, Matrix& res){
        const long long mod = x.mod, limit = mod * mod;  /// sums stay below 2 * mod^2 < 2^63

        for (int i = 0; i < x.n; i++){
            vector<long long>& acc = res.a[i];
            fill(acc.begin(), acc.end(), 0);
            for (int k = 0; k < x.m; k++){
                long long v = x.a[i][k] % mod;
                if (v < 0) v += mod;
                if (!v) continue;
                const vector<long long>& row = y.a[k];
                for (int j = 0; j < y.m; j++){
                    acc[j] += v * row[j];
                    if (acc[j] >= limit) acc[j] -= limit;
                }
            }
            for (long long& value : acc) value %= mod;
        }
    }
};

const long long MOD = 1000000000;

int k, n;
long long b[12], c[12];

long long solve(){
    if (n <= k) return b[n] % MOD;

    Matrix M(k, k, MOD);
    for (int j = 0; j < k; j++) M[0][j] = c[j + 1];
    for (int i = 1; i < k; i++) M[i][i - 1] = 1;

    Matrix P = M.pow(n - k);
    long long res = 0;
    for (int j = 0; j < k; j++) res = (res + P[0][j] * (b[k - j] % MOD)) % MOD;
    return res;
}

int main(){
    int t, i;

    scanf("%d", &t);
    while (t--){
        scanf("%d", &k);
        for (i = 1; i <= k; i++) scanf("%lld", &b[i]);
        for (i = 1; i <= k; i++) scanf("%lld", &c[i]);
        scanf("%d", &n);
        printf("%lld\n", solve());
    }

    return 0;
}
