/***
 *
 * Matrix Permanent
 * Permanent of an integer matrix modulo any m by Ryser's formula, and the exact permanent when it is not too large
 *
 * Complexity: O(2^n * n), O(n^2) extra space
 *
 * perm(a) = sum over permutations p of a[0][p[0]] * a[1][p[1]] * ... * a[n - 1][p[n - 1]]
 * Ryser: perm(a) = sum over column subsets S of (-1)^(n - |S|) * prod over rows i of (sum of a[i][j] for j in S)
 * Subsets are visited in Gray code order, so each step adds or removes one column from the row sums in O(n)
 *
 * permanent_mod(a, m): perm(a) mod m in [0, m) for any modulus 1 <= m < 2^62, prime or composite
 * permanent(a): the exact signed permanent, correct whenever |perm(a)| <= 2258704744122758558 (~2.26e18)
 *     computed modulo 4517409488245517117 and mapped back to the nearest signed value
 *
 * Entries may be any long long, negative included, a 0 x 0 matrix has permanent 1
 * n = 20 takes about 2 * 10^7 modular multiplications
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long permanent_mod(const vector<vector<long long>>& a, long long m){
    int n = a.size();
    assert(1 <= m && m < (1LL << 62) && n < 63);
    vector<vector<long long>> b(n, vector<long long>(n));
    for (int i = 0; i < n; i++){
        assert((int)a[i].size() == n);
        for (int j = 0; j < n; j++) b[i][j] = (a[i][j] % m + m) % m;
    }

    vector<long long> row_sum(n, 0);
    long long res = n == 0 ? 1 % m : 0;
    for (long long k = 1; k < (1LL << n); k++){
        long long gray = k ^ (k >> 1);
        int col = __builtin_ctzll(k);
        bool added = gray >> col & 1;
        for (int i = 0; i < n; i++){
            row_sum[i] += added ? b[i][col] : m - b[i][col];
            if (row_sum[i] >= m) row_sum[i] -= m;
        }

        long long prod = 1;
        for (int i = 0; i < n && prod; i++) prod = (long long)((__int128)prod * row_sum[i] % m);
        if ((n - __builtin_popcountll(gray)) % 2) prod = prod ? m - prod : 0;
        res += prod;
        if (res >= m) res -= m;
    }

    return res;
}

long long permanent(const vector<vector<long long>>& a){
    const long long P = 4517409488245517117LL;
    long long r = permanent_mod(a, P);
    return r > P / 2 ? r - P : r;
}

int main(){
    assert(permanent({{1, 2}, {3, 4}}) == 10);  /// 1 * 4 + 2 * 3, a determinant without the signs
    assert(permanent({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}) == 450);
    assert(permanent({{-1, 2}, {3, -4}}) == 10);
    assert(permanent({}) == 1);

    /// J - I counts derangements: 9 of 4 elements
    assert(permanent({{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}}) == 9);

    assert(permanent_mod({{1, 2}, {3, 4}}, 7) == 3);
    vector<vector<long long>> ones(20, vector<long long>(20, 1));
    assert(permanent_mod(ones, 1000000007) == 146326063);  /// 20! mod 1e9 + 7
    return 0;
}
