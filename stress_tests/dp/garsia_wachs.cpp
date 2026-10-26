#include "../common.h"

#define main library_main
#include "../../code_library/dp/garsia_wachs.cpp"
#undef main

#define main knuth_main
#include "../../code_library/dp/knuth_optimization.cpp"
#undef main

long long naive_merge_cost(const vector<long long>& weights){
    int n = weights.size();
    if (n == 0) return 0;
    vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + weights[i];

    vector<vector<long long>> dp(n, vector<long long>(n, 0));
    for (int len = 2; len <= n; len++){
        for (int i = 0, j = len - 1; j < n; i++, j++){
            dp[i][j] = LLONG_MAX;
            for (int k = i; k < j; k++) dp[i][j] = min(dp[i][j], dp[i][k] + dp[k + 1][j]);
            dp[i][j] += pre[j + 1] - pre[i];
        }
    }
    return dp[0][n - 1];
}

long long knuth_merge_cost(const vector<long long>& weights){
    int n = weights.size();
    vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + weights[i];
    return knuth_dp<long long>(n, [&](int i, int j){ return pre[j] - pre[i]; })[0][n];
}

vector<long long> shaped(int n, int shape, long long max_weight){
    vector<long long> w(n);
    for (int i = 0; i < n; i++){
        if (shape <= 2) w[i] = stress::rand_int(0, max_weight);
        if (shape == 3) w[i] = max_weight;
        if (shape == 4) w[i] = stress::rand_int(0, 1) ? 0 : stress::rand_int(0, max_weight);
        if (shape == 5) w[i] = min(i, n - 1 - i);  /// a peak, merged piles cross long runs of smaller piles
        if (shape == 6) w[i] = abs(n / 2 - i);  /// a valley
    }
    if (shape == 1) sort(w.begin(), w.end());
    if (shape == 2) sort(w.rbegin(), w.rend());
    return w;
}

/// Exhaustive weights 0..3 for n <= 8 hit every tie and zero pattern the <= and < comparisons decide
int main(){
    for (int n = 0; n <= 8; n++){
        int combos = 1;
        for (int i = 0; i < n; i++) combos *= 4;
        for (int mask = 0; mask < combos; mask++){
            vector<long long> w(n);
            for (int i = 0, m = mask; i < n; i++, m /= 4) w[i] = m % 4;
            assert(garsia_wachs(w) == naive_merge_cost(w));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 60);
        long long max_weight = stress::rand_int(0, 1) ? 5 : 1000000000;
        vector<long long> w = shaped(n, stress::rand_int(0, 6), max_weight);
        assert(garsia_wachs(w) == naive_merge_cost(w));
    }

    for (long long it = 0; it < stress::scaled(40); it++){
        int n = stress::rand_int(61, 1500);
        long long max_weight = stress::rand_int(0, 1) ? 3 : 1000000000000LL;
        vector<long long> w = shaped(n, it % 7, max_weight);
        assert(garsia_wachs(w) == knuth_merge_cost(w));
    }

    /// The header's range, total * ceil(log2 n) < LLONG_MAX: equal weights on 2^k piles build the full tree
    for (int k = 1; k <= 17; k++){
        int n = 1 << k;
        long long w = LLONG_MAX / k / n;
        assert(garsia_wachs(vector<long long>(n, w)) == w * n * k);
    }
    assert(garsia_wachs(vector<long long>{LLONG_MAX / 2, LLONG_MAX / 2}) == LLONG_MAX - 1);
    assert(garsia_wachs(vector<long long>{LLONG_MAX - 1}) == 0);
    vector<long long> heavy = shaped(1000, 0, LLONG_MAX / 10 / 1000);
    assert(garsia_wachs(heavy) == knuth_merge_cost(heavy));

    assert(garsia_wachs(vector<int>{3, 1, 4, 1, 5}) == naive_merge_cost({3, 1, 4, 1, 5}));

    return 0;
}
