/***
 *
 * Matrix
 * Matrix multiplication and exponentiation, modulo mod or over any semiring, plus transpose and rotation
 *
 * Complexity:
 *   - multiplication: O(n m p) for an n x m times an m x p matrix
 *   - pow(e): O(n^3 log e)
 *
 * Matrix a(n, m, mod): n x m zero matrix of long long modulo mod, 1 <= mod < 2^31
 * Matrix<T, Ops = PlusTimes<T>> a(n, m): n x m matrix of Ops::zero(), no modulus, like Matrix<long double> a(n, m)
 * Matrix<T, Ops>::identity(n, mod = 0): the n x n identity, Matrix<>::identity(n, mod) for the modular one
 * a[i][j] = value: any T, with a modulus any long long, reduced modulo mod when multiplying
 * a * b: the product, needs a.m == b.n and the same mod
 * a.is_normalized(): whether every entry is already in [0, mod), always true without a modulus
 * a.normalized(): a copy with every entry in [0, mod), a plain copy without a modulus
 * a.pow(e): a^e for a square matrix and e >= 0, a.pow(0) is the identity
 * a.rotate(): the matrix turned 90 degrees clockwise
 * a.transpose(): the m x n transpose
 *
 * Notes:
 *   - Recurrence f(n) = 2 f(n - 1) + 3 f(n - 2): with M = {{2, 3}, {1, 0}}, (f(n + 1), f(n)) = M^n * (f(1), f(0))
 *   - Ops holds static zero(), one(), add(x, y), mul(x, y) of a semiring: zero absorbs in mul, mul distributes over add
 *   - Min-plus: zero = INF, one = 0, add = min, mul = + capped at INF, then a.pow(k)[u][v] is the shortest k-edge walk
 *   - A modulus needs the default Matrix<long long>, so a custom Ops can never be silently replaced by modular math
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
    static constexpr bool MODULAR = is_same_v<T, long long> && is_same_v<Ops, PlusTimes<T>>;

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

int main(){
    /// Fibonacci by exponentiation, then a plain product
    const long long MOD = 1000000007;

    Matrix fib(2, 2, MOD);
    fib.a = {{1, 1}, {1, 0}};
    assert(fib.pow(10)[0][1] == 55);  /// F(10)

    Matrix a(2, 3, MOD), b(3, 2, MOD);
    a.a = {{1, 2, 3}, {4, 5, 6}};
    b.a = {{7, 8}, {9, 10}, {11, 12}};
    assert(((a * b).a == vector<vector<long long>>{{58, 64}, {139, 154}}));

    /// Negative entries are reduced
    Matrix neg(1, 1, 7);
    neg[0][0] = -1;
    assert(!neg.is_normalized());
    assert(neg.normalized()[0][0] == 6);
    assert(neg.pow(3)[0][0] == 6);  /// (-1)^3 = -1 = 6 mod 7

    assert((a.transpose().a == vector<vector<long long>>{{1, 4}, {2, 5}, {3, 6}}));
    assert((a.rotate().a == vector<vector<long long>>{{4, 1}, {5, 2}, {6, 3}}));  /// clockwise

    /// Real entries: halving twice
    Matrix<long double> half(1, 1);
    half[0][0] = 0.5;
    assert(half.pow(2)[0][0] == 0.25);

    /// Min-plus: shortest walk with exactly 2 edges
    struct MinPlus{
        static long long zero(){ return LLONG_MAX; }
        static long long one(){ return 0; }
        static long long add(long long x, long long y){ return min(x, y); }
        static long long mul(long long x, long long y){ return x == LLONG_MAX || y == LLONG_MAX ? LLONG_MAX : x + y; }
    };

    Matrix<long long, MinPlus> dist(3, 3);
    dist[0][1] = 4, dist[1][2] = 1, dist[0][2] = 9;
    assert(dist.pow(2)[0][2] == 5);  /// 0 -> 1 -> 2
    assert(dist.pow(2)[0][1] == LLONG_MAX);  /// no 2-edge walk from 0 to 1

    return 0;
}
