#include "common.h"

#define main library_main
#include "../code_library/determinant.cpp"
#undef main

/// Exact determinant by Laplace expansion along the first row, in __int128
__int128 laplace(const vector<vector<long long>>& a){
    int n = a.size();
    if (n == 0) return 1;
    __int128 res = 0;
    for (int c = 0; c < n; c++){
        vector<vector<long long>> minor;
        for (int r = 1; r < n; r++){
            vector<long long> row;
            for (int k = 0; k < n; k++){
                if (k != c) row.push_back(a[r][k]);
            }
            minor.push_back(row);
        }
        __int128 sub = laplace(minor) * a[0][c];
        res += c % 2 ? -sub : sub;
    }
    return res;
}

/// Exact determinant by fraction free Bareiss elimination, fine while the minors stay small
__int128 bareiss(vector<vector<__int128>> a){
    int n = a.size(), sign = 1;
    __int128 prev = 1;
    for (int i = 0; i < n; i++){
        if (a[i][i] == 0){
            int p = i + 1;
            while (p < n && a[p][i] == 0) p++;
            if (p == n) return 0;
            swap(a[i], a[p]), sign = -sign;
        }
        for (int j = i + 1; j < n; j++){
            for (int k = i + 1; k < n; k++) a[j][k] = (a[j][k] * a[i][i] - a[j][i] * a[i][k]) / prev;
        }
        prev = a[i][i];
    }
    return sign * (n ? a[n - 1][n - 1] : 1);
}

long long reduce(__int128 x, long long m){
    x %= m;
    return (long long)(x < 0 ? x + m : x);
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int n = stress::rand_int(0, 5);
        long long range = it % 3 == 0 ? 3 : (it % 3 == 1 ? 1000 : 1000000);
        vector<vector<long long>> a(n, vector<long long>(n));
        for (auto& row : a){
            for (auto& x : row) x = stress::rand_int(-range, range);
        }
        if (n >= 2 && it % 7 == 0) a[1] = a[0];
        __int128 exact = laplace(a);

        long long m = it % 4 == 0 ? stress::rand_int(1, 30) : (it % 4 == 1 ? (1LL << 62) - 1 : stress::rand_int(1, (1LL << 62) - 1));
        assert(determinant_mod(a, m) == reduce(exact, m));
        if (exact >= -2258704744122758558LL && exact <= 2258704744122758558LL) assert(determinant(a) == (long long)exact);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(6, 18);
        vector<vector<long long>> a(n, vector<long long>(n));
        vector<vector<__int128>> b(n, vector<__int128>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) a[i][j] = b[i][j] = stress::rand_int(-2, 2);
        }
        __int128 exact = bareiss(b);
        assert(determinant(a) == (long long)exact);
        long long m = stress::rand_int(1, 1000000007);
        assert(determinant_mod(a, m) == reduce(exact, m));
    }

    const long long HALF = 2258704744122758558LL;
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long v = it < 6 ? vector<long long>{HALF, -HALF, HALF - 1, 1 - HALF, 0, 1}[it] : stress::rand_int(-HALF, HALF);
        assert(determinant({{v}}) == v);
        assert(determinant({{v, 0}, {0, 1}}) == v && determinant({{0, v}, {1, 0}}) == -v);
    }
    return 0;
}
