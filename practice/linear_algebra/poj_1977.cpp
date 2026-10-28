/***
 *
 * Matrix
 * Modular matrix multiplication and exponentiation, plus transpose and rotation
 *
 * Complexity:
 *   - multiplication: O(n m p) for an n x m times an m x p matrix
 *   - pow(e): O(n^3 log e), with no allocation inside the loop
 *
 * Matrix a(n, m, mod): n x m zero matrix modulo mod, 1 <= mod < 2^31
 * Matrix::identity(n, mod): the n x n identity modulo mod
 * a[i][j] = value: any long long, reduced modulo mod when multiplying
 * a * b: the product, needs a.m == b.n and the same mod
 * a.pow(e): a^e for a square matrix and e >= 0, a.pow(0) is the identity
 * a.normalized(): a copy with every entry in [0, mod)
 * a.is_normalized(): whether every entry is already in [0, mod)
 * a.transpose(): the m x n transpose
 * a.rotate(): the matrix turned 90 degrees clockwise
 *
 * Linear recurrence f(n) = 2 f(n - 1) + 3 f(n - 2): with M = {{2, 3}, {1, 0}}, (f(n + 1), f(n)) = M^n * (f(1), f(0))
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

    bool is_normalized() const{
        for (const auto& row : a){
            for (long long x : row){
                if (x < 0 || x >= mod) return false;
            }
        }

        return true;
    }

    Matrix operator*(const Matrix& other) const{
        assert(m == other.n && mod == other.mod);

        Matrix res(n, other.m, mod);
        if (other.is_normalized()) multiply(*this, other, res);
        else multiply(*this, other.normalized(), res);
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

    /// res = x * y in res's existing rows; y must be normalized and res must not alias x or y
    static void multiply(const Matrix& x, const Matrix& y, Matrix& res){
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

int index_of(map<string, int>& id, const string& name){
    map<string, int>::iterator it = id.find(name);
    if (it != id.end()) return it->second;
    int k = id.size();
    id[name] = k;
    return k;
}

int main(){
    int tests;

    cin >> tests;
    while (tests--){
        int n;
        long long t;
        cin >> n >> t;

        map<string, int> id;
        Matrix step = Matrix::identity(n, 2), marks(n, 1, 2);
        for (int i = 0; i < n; i++){
            string name;
            int chalk, k;
            cin >> name >> chalk >> k;
            int u = index_of(id, name);
            marks[u][0] = chalk;
            while (k--){
                cin >> name;
                step[index_of(id, name)][u]++;
            }
        }

        Matrix last = step.pow(t - 1) * marks;
        int winners = 0;
        for (int i = 0; i < n; i++) winners += last[i][0];
        cout << winners << endl;
    }

    return 0;
}
