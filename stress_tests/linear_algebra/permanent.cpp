#include "../common.h"

#define main library_main
#include "../../code_library/linear_algebra/permanent.cpp"
#undef main

/// Exact permanent as the sum over all permutations, in __int128
__int128 brute_exact(const vector<vector<long long>>& a){
    int n = a.size();
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    __int128 res = 0;

    do{
        __int128 prod = 1;
        for (int i = 0; i < n; i++) prod *= a[i][p[i]];
        res += prod;
    } while (next_permutation(p.begin(), p.end()));

    return res;
}

long long reduce(__int128 x, long long m){
    x %= m;
    return (long long)(x < 0 ? x + m : x);
}

/// Permanent mod m as the sum over all permutations
long long brute_mod(const vector<vector<long long>>& a, long long m){
    int n = a.size();
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    vector<vector<long long>> b = a;
    for (auto& row : b){
        for (auto& x : row) x = reduce(x, m);
    }
    __int128 res = 0;

    do{
        __int128 prod = 1 % m;
        for (int i = 0; i < n; i++) prod = prod * b[i][p[i]] % m;
        res = (res + prod) % m;
    } while (next_permutation(p.begin(), p.end()));

    return (long long)res;
}

/// Exact permanent by the row-by-row subset DP of the original source, dp[mask] = ways to finish rows popcount(mask)..n-1
__int128 subset_dp(const vector<vector<long long>>& a){
    int n = a.size();
    vector<__int128> dp(1 << n, 0);
    dp[(1 << n) - 1] = 1;

    for (int mask = (1 << n) - 2; mask >= 0; mask--){
        int i = __builtin_popcount(mask);
        for (int j = 0; j < n; j++){
            if (!(mask >> j & 1)) dp[mask] += dp[mask | 1 << j] * a[i][j];
        }
    }

    return dp[0];
}

long long random_modulus(long long it){
    if (it % 4 == 0) return stress::rand_int(1, 30);
    if (it % 4 == 1) return (1LL << 62) - 1;
    return stress::rand_int(1, (1LL << 62) - 1);
}

int main(){
    const long long HALF = 2258704744122758558LL;

    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = it % 16 == 15 ? 8 : stress::rand_int(0, 7);
        long long range = it % 3 == 0 ? 3 : (it % 3 == 1 ? 1000 : 1000000000000000000LL);
        vector<vector<long long>> a(n, vector<long long>(n));
        for (auto& row : a){
            for (auto& x : row) x = stress::rand_int(-range, range);
        }
        if (n >= 2 && it % 7 == 0) a[1] = a[0];
        if (n >= 1 && it % 11 == 0) a[0].assign(n, 0);

        long long m = random_modulus(it);
        assert(permanent_mod(a, m) == brute_mod(a, m));
        if (range <= 1000){
            __int128 exact = brute_exact(a);
            if (exact >= -HALF && exact <= HALF) assert(permanent(a) == (long long)exact);
        }
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(9, it % 10 == 0 ? 16 : 13);
        vector<vector<long long>> a(n, vector<long long>(n));
        for (auto& row : a){
            for (auto& x : row) x = stress::rand_int(-2, 2);
        }

        __int128 exact = subset_dp(a);
        assert(permanent(a) == (long long)exact);
        long long m = random_modulus(it);
        assert(permanent_mod(a, m) == reduce(exact, m));
    }

    for (long long it = 0; it < stress::scaled(2000); it++){
        long long v = it < 6 ? vector<long long>{HALF, -HALF, HALF - 1, 1 - HALF, 0, 1}[it] : stress::rand_int(-HALF, HALF);
        assert(permanent({{v}}) == v);
        assert(permanent({{v, 0}, {0, 1}}) == v && permanent({{0, v}, {-1, 0}}) == -v);
    }

    return 0;
}
