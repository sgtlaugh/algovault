/***
 *
 * Characteristic Polynomial modulo a prime
 * Coefficients of det(xI - A) for a square matrix A over Z/pZ
 *
 * Complexity: O(n^3) time, O(n^2) memory
 *
 * characteristic_polynomial(a, p) returns c[0..n] with det(xI - a) = c[0] + c[1] x + ... + c[n] x^n, c[n] = 1
 *     p must be a prime, 2 <= p < 2^31, entries may be any long long and are reduced into [0, p)
 *     a 0 x 0 matrix gives {1}
 *
 * A similarity transform first reduces a to upper Hessenberg form (zero below the subdiagonal),
 * then the polynomials of its leading principal submatrices follow by expanding along the last column
 * No evaluation points or interpolation are used, so any p works, including p <= n
 *
 * det(a) = (-1)^n c[0] and trace(a) = -c[n - 1]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long expo_mod(long long a, long long e, long long p){
    long long res = 1;

    while (e){
        if (e & 1) res = res * a % p;
        a = a * a % p;
        e >>= 1;
    }

    return res;
}

vector<long long> characteristic_polynomial(vector<vector<long long>> a, long long p){
    int n = a.size();
    assert(2 <= p && p < (1LL << 31));
    for (auto& row : a){
        assert((int)row.size() == n);
        for (auto& x : row) x = (x % p + p) % p;
    }

    /// sums of products are kept below p^2 < 2^62 by subtraction, a % per term doubled the runtime at n = 500
    const long long psq = p * p;
    for (int c = 0; c + 2 < n; c++){
        int pivot = c + 1;
        while (pivot < n && a[pivot][c] == 0) pivot++;
        if (pivot == n) continue;

        /// every row operation is paired with the inverse column operation, so the matrix stays similar to a
        swap(a[pivot], a[c + 1]);
        for (int i = 0; i < n; i++) swap(a[i][pivot], a[i][c + 1]);
        long long inv = expo_mod(a[c + 1][c], p - 2, p);
        vector<long long> f(n, 0);
        for (int r = c + 2; r < n; r++){
            f[r] = a[r][c] * inv % p;
            if (f[r] == 0) continue;
            for (int j = c; j < n; j++) a[r][j] = (a[r][j] + (p - f[r]) * a[c + 1][j]) % p;
        }

        /// the row operations of one column commute, so their column operations batch into one row-major pass
        for (int i = 0; i < n; i++){
            long long sum = a[i][c + 1];
            for (int r = c + 2; r < n; r++){
                sum += f[r] * a[i][r];
                if (sum >= psq) sum -= psq;
            }
            a[i][c + 1] = sum % p;
        }
    }

    vector<vector<long long>> poly(n + 1);
    poly[0] = {1};
    for (int k = 0; k < n; k++){
        vector<long long>& cur = poly[k + 1];
        cur.assign(k + 2, 0);
        for (int j = 0; j <= k; j++){
            cur[j + 1] = (cur[j + 1] + poly[k][j]) % p;
            cur[j] = (cur[j] + (p - a[k][k]) * poly[k][j]) % p;
        }

        long long t = 1;
        for (int i = k - 1; i >= 0 && t != 0; i--){
            t = t * a[i + 1][i] % p;
            long long f = t * a[i][k] % p;
            for (int j = 0; j <= i; j++) cur[j] = (cur[j] + (p - f) * poly[i][j]) % p;
        }
    }

    return poly[n];
}

int main(){
    const long long mod = 998244353;
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& v : row) scanf("%lld", &v);
    }

    vector<long long> c = characteristic_polynomial(a, mod);
    for (int i = 0; i <= n; i++) printf("%lld%c", c[i], i == n ? '\n' : ' ');
    return 0;
}
