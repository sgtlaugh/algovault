/***
 *
 * Integer Determinant
 * Determinant of an integer matrix modulo any m, and the exact determinant when it is not too large
 *
 * Complexity: O(n^3 + n^2 log m)
 *
 * determinant_mod(a, m): det(a) mod m in [0, m) for any modulus 1 <= m < 2^62, prime or composite
 *     row reduction by repeated division, like the Euclidean algorithm, so no modular inverses are needed
 * determinant(a): the exact signed determinant, correct whenever |det(a)| <= 2258704744122758558 (~2.26e18)
 *     computed modulo the prime 4517409488245517117 and mapped back to the nearest signed value
 *
 * Entries may be any long long, negative included, a 0 x 0 matrix has determinant 1
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

long long determinant_mod(vector<vector<long long>> a, long long m){
    int n = a.size();
    assert(1 <= m && m < (1LL << 62));
    for (auto& row : a){
        assert((int)row.size() == n);
        for (auto& x : row) x %= m;
    }

    __int128 res = 1;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            while (a[j][i] != 0){
                long long t = a[i][i] / a[j][i];
                if (t){
                    for (int k = i; k < n; k++) a[i][k] = (long long)((a[i][k] - (__int128)a[j][k] * t) % m);
                }
                swap(a[i], a[j]);
                res = -res;
            }
        }
        res = res * a[i][i] % m;
        if (res == 0) return 0;
    }

    res %= m;
    return (long long)(res < 0 ? res + m : res);
}

long long determinant(const vector<vector<long long>>& a){
    const long long P = 4517409488245517117LL;
    long long r = determinant_mod(a, P);
    return r > P / 2 ? r - P : r;
}

int main(){
    assert(determinant({{1, 2}, {3, 4}}) == -2);
    assert(determinant({{2, 0, 0}, {0, 3, 0}, {0, 0, 5}}) == 30);
    assert(determinant({{0, 1}, {1, 0}}) == -1);
    assert(determinant({{1, 2}, {2, 4}}) == 0);
    assert(determinant({{6, 1, 1}, {4, -2, 5}, {2, 8, 7}}) == -306);
    assert(determinant({}) == 1);
    assert(determinant({{-7}}) == -7);

    assert(determinant_mod({{1, 2}, {3, 4}}, 7) == 5);
    assert(determinant_mod({{2, 1}, {1, 3}}, 6) == 5);
    assert(determinant_mod({{2, 0}, {0, 3}}, 6) == 0);
    assert(determinant_mod({{4, 1}, {1, 4}}, 1) == 0);
    assert(determinant_mod({}, 10) == 1);

    const long long HALF = 2258704744122758558LL;
    assert(determinant({{HALF}}) == HALF);
    assert(determinant({{-HALF}}) == -HALF);
    assert(determinant({{1000000000, 0}, {0, 1000000000}}) == 1000000000000000000LL);

    return 0;
}
