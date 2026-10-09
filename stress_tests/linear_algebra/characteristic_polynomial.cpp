#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/characteristic_polynomial.cpp"
#undef main

long long inv_mod(long long a, long long p){
    long long res = 1, e = p - 2;

    while (e){
        if (e & 1) res = res * a % p;
        a = a * a % p;
        e >>= 1;
    }

    return res;
}

/// Determinant mod p by Gaussian elimination with modular inverses
long long det_mod(vector<vector<long long>> a, long long p){
    int n = a.size();
    long long res = 1;

    for (int c = 0; c < n; c++){
        int r = c;
        while (r < n && a[r][c] == 0) r++;
        if (r == n) return 0;
        if (r != c) swap(a[r], a[c]), res = p - res;
        res = res * a[c][c] % p;

        long long inv = inv_mod(a[c][c], p);
        for (int i = c + 1; i < n; i++){
            long long f = a[i][c] * inv % p;
            for (int j = c; j < n; j++) a[i][j] = (a[i][j] - f * a[c][j] % p + p) % p;
        }
    }

    return res % p;
}

/// det(xI - a) at x = 0..n, then Lagrange interpolation, needs p > n
vector<long long> by_interpolation(const vector<vector<long long>>& a, long long p){
    int n = a.size();
    vector<long long> ys(n + 1);
    for (int x = 0; x <= n; x++){
        vector<vector<long long>> b(n, vector<long long>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) b[i][j] = ((i == j ? x : 0) - a[i][j] % p + 2 * p) % p;
        }
        ys[x] = det_mod(b, p);
    }

    vector<long long> res(n + 1, 0);
    for (int i = 0; i <= n; i++){
        vector<long long> num = {1};
        long long den = 1;
        for (int j = 0; j <= n; j++){
            if (j == i) continue;
            vector<long long> next(num.size() + 1, 0);
            for (int k = 0; k < (int)num.size(); k++){
                next[k + 1] = (next[k + 1] + num[k]) % p;
                next[k] = (next[k] + (p - j) * num[k]) % p;
            }
            num = next;
            den = den * ((i - j + p) % p) % p;
        }
        long long scale = ys[i] * inv_mod(den, p) % p;
        for (int k = 0; k <= n; k++) res[k] = (res[k] + num[k] * scale) % p;
    }

    return res;
}

/// Leibniz expansion of det(xI - a) over all permutations with polynomial entries, any p
vector<long long> by_permutations(const vector<vector<long long>>& a, long long p){
    int n = a.size();
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    vector<long long> res(n + 1, 0);

    do{
        int inversions = 0;
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++) inversions += perm[i] > perm[j];
        }

        vector<long long> prod = {1};
        for (int i = 0; i < n; i++){
            long long c = (p - a[i][perm[i]] % p) % p;
            vector<long long> next(prod.size() + 1, 0);
            for (int k = 0; k < (int)prod.size(); k++){
                next[k] = (next[k] + c * prod[k]) % p;
                if (perm[i] == i) next[k + 1] = (next[k + 1] + prod[k]) % p;
            }
            prod = next;
        }
        for (int k = 0; k <= n; k++){
            long long v = prod[k] % p;
            res[k] = (res[k] + (inversions % 2 ? p - v : v)) % p;
        }
    } while (next_permutation(perm.begin(), perm.end()));

    return res;
}

vector<vector<long long>> random_matrix(int n, long long p){
    long long range = stress::rand_int(0, 2) == 0 ? 2 : p - 1;
    int zero_pct = stress::rand_int(0, 2) == 0 ? 80 : (stress::rand_int(0, 1) ? 40 : 0);
    vector<vector<long long>> a(n, vector<long long>(n));
    for (auto& row : a){
        for (auto& x : row) x = stress::rand_int(0, 99) < zero_pct ? 0 : stress::rand_int(-range, range);
    }
    if (n >= 2 && stress::rand_int(0, 6) == 0) a[1] = a[0];
    return a;
}

int main(){
    const vector<long long> small_primes = {2, 3, 5, 7};
    const vector<long long> big_primes = {998244353, 1000000007, 2147483647};

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 6);
        long long p = it % 2 ? small_primes[stress::rand_int(0, 3)] : big_primes[stress::rand_int(0, 2)];
        auto a = random_matrix(n, p);
        assert(characteristic_polynomial(a, p) == by_permutations(a, p));
    }

    for (long long it = 0; it < stress::scaled(200); it++){
        int n = it < 30 ? it : stress::rand_int(1, 45);
        long long p = it % 3 == 0 ? 2147483647 : (it % 3 == 1 ? 998244353 : 1000003);
        auto a = random_matrix(n, p);
        if (it % 5 == 0){
            for (auto& row : a){
                for (auto& x : row) x = stress::rand_int(LLONG_MIN / 2, LLONG_MAX / 2);
            }
        }
        vector<vector<long long>> reduced = a;
        for (auto& row : reduced){
            for (auto& x : row) x = (x % p + p) % p;
        }
        assert(characteristic_polynomial(a, p) == by_interpolation(reduced, p));
    }

    return 0;
}
