#include "../common.h"

#define main library_main
#include "../../code_library/dp/concave_1d1d_dp.cpp"
#undef main

/// O(n^2) DP over every i < x, on general random Monge costs and convex-gap costs
template <typename Cost>
vector<long long> brute(int n, long long dp0, Cost w){
    vector<long long> dp(n + 1, dp0);
    for (int x = 1; x <= n; x++){
        dp[x] = LLONG_MAX;
        for (int i = 0; i < x; i++) dp[x] = min(dp[x], dp[i] + w(i, x));
    }

    return dp;
}

/// w(i, x) = r[i] + c[x] + sum of D[p][q] >= 0 over p <= i, q >= x satisfies the quadrangle inequality,
/// and every Monge matrix has this form, so this samples the whole valid input space
vector<vector<long long>> random_monge(int n, int density, int range){
    vector<vector<long long>> D(n + 2, vector<long long>(n + 2, 0)), S(n + 2, vector<long long>(n + 2, 0));
    for (int p = 0; p <= n; p++){
        for (int q = 0; q <= n; q++){
            if (stress::rand_int(1, 100) <= density) D[p][q] = stress::rand_int(0, range);
        }
    }

    for (int p = 0; p <= n; p++){
        for (int q = n; q >= 0; q--) S[p][q] = D[p][q] + (p ? S[p - 1][q] : 0) + S[p][q + 1] - (p ? S[p - 1][q + 1] : 0);
    }

    vector<long long> r(n + 1), c(n + 1);
    for (auto& v : r) v = stress::rand_int(-range * 3, range * 3);
    for (auto& v : c) v = stress::rand_int(-range * 3, range * 3);
    for (int i = 0; i <= n; i++){
        for (int x = 0; x <= n; x++) S[i][x] += r[i] + c[x];
    }

    return S;
}

/// g(0..n) with nondecreasing differences, so g(x - i) is a convex-gap cost
vector<long long> random_convex(int n, int range){
    vector<long long> g(n + 1);
    long long slope = stress::rand_int(-range * n, range), value = stress::rand_int(-range, range);
    for (int d = 0; d <= n; d++){
        g[d] = value;
        slope += stress::rand_int(0, range);
        value += slope;
    }

    return g;
}

template <typename Cost>
void check(int n, long long dp0, Cost w){
    vector<int> opt;
    auto dp = concave_1d1d_dp(n, dp0, w, &opt);
    assert(dp == brute(n, dp0, w));
    assert(opt[0] == -1);
    for (int x = 1; x <= n; x++) assert(0 <= opt[x] && opt[x] < x && dp[opt[x]] + w(opt[x], x) == dp[x]);
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = it < 40 ? it % 20 : stress::rand_int(0, 40);
        int density = vector<int>{2, 10, 50, 100}[it % 4], range = it % 3 ? 3 : 1000;
        auto M = random_monge(n, density, range);
        check(n, stress::rand_int(-50, 50), [&](int i, int x){ return M[i][x]; });
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, it < 280 ? 300 : 2000), range = it % 2 ? 2 : 100000;
        auto g = random_convex(n, range);
        vector<long long> A(n + 1), B(n + 1);
        for (int i = 1; i <= n; i++) A[i] = A[i - 1] + stress::rand_int(0, 50), B[i] = B[i - 1] + stress::rand_int(0, 50);
        long long scale = stress::rand_int(0, 1);
        check(n, 0, [&](int i, int x){ return g[x - i] - scale * A[x] * B[i] + B[x] * 7 - A[i]; });
    }

    int n = 100000;
    long long calls = 0;
    auto dp = concave_1d1d_dp(n, 0LL, [&](int i, int x){ calls++; return 1LL * (x - i) * (x - i) + 10000; });
    assert(calls <= 4LL * n * 17);
    assert(dp[n] == 200LL * n);

    return 0;
}
