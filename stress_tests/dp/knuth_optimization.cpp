#include "../common.h"

#define main library_main
#include "../../code_library/dp/knuth_optimization.cpp"
#undef main

template <typename Cost>
vector<vector<long long>> naive_dp(int n, Cost cost){
    vector<vector<long long>> dp(n + 1, vector<long long>(n + 1, 0));
    for (int len = 2; len <= n; len++){
        for (int i = 0, j = len; j <= n; i++, j++){
            dp[i][j] = LLONG_MAX;
            for (int k = i + 1; k < j; k++) dp[i][j] = min(dp[i][j], dp[i][k] + dp[k][j]);
            dp[i][j] += cost(i, j);
        }
    }
    return dp;
}

/// Shapes 0-3 are f(interval sum of nonnegative weights) with f convex and nondecreasing, the class the header claims
/// Shape 4 is the optimal BST cost, shape 5 adds a squared length term that no interval sum produces
long long shaped_cost(int shape, const vector<long long>& pre, int i, int j, long long threshold){
    long long s = pre[j] - pre[i];
    if (shape == 0) return s;
    if (shape == 1) return s * s;
    if (shape == 2) return max(0LL, s - threshold);
    if (shape == 3) return s * s + 3 * s + threshold;
    if (shape == 4) return pre[j - 1] - pre[i];
    return (long long)(j - i) * (j - i) + s;
}

void check(const vector<long long>& weights, int shape, long long threshold){
    int n = weights.size();
    vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + weights[i];

    auto cost = [&](int i, int j){ return shaped_cost(shape, pre, i, j, threshold); };
    assert((knuth_dp<long long>(n, cost) == naive_dp(n, cost)));
}

/// Knuth DP tables against the O(n^3) recurrence over every interval, exhaustive weights for small n, random beyond
int main(){
    for (int n = 0; n <= 7; n++){
        int combos = 1;
        for (int i = 0; i < n; i++) combos *= 3;
        for (int mask = 0; mask < combos; mask++){
            vector<long long> weights(n);
            for (int i = 0, m = mask; i < n; i++, m /= 3) weights[i] = m % 3;
            for (int shape = 0; shape < 6; shape++) check(weights, shape, 2);
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(0, it < 1400 ? 40 : 200);
        long long max_weight = stress::rand_int(0, 1) ? 3 : 100000;
        vector<long long> weights(n);
        for (auto& w : weights) w = stress::rand_int(0, max_weight);
        check(weights, stress::rand_int(0, 5), stress::rand_int(0, max_weight * 4));
    }

    return 0;
}
