#include "../common.h"

#define main library_main
#include "../../code_library/dp/aliens_trick.cpp"
#undef main

const long long INF = LLONG_MAX / 4;

/// Convex f on [first, n] given directly; the solver scans every count, so ties come from repeated slopes
void check_synthetic(){
    int first = stress::rand_int(0, 3), n = first + stress::rand_int(0, 10);
    int range = stress::rand_int(1, 2) == 1 ? 3 : 1000;

    vector<long long> slopes;
    for (int c = first + 1; c <= n; c++) slopes.push_back(stress::rand_int(-range, range));
    sort(slopes.begin(), slopes.end());

    vector<long long> f(n + 1, INF);
    f[first] = stress::rand_int(-range, range);
    for (int c = first + 1; c <= n; c++) f[c] = f[c - 1] + slopes[c - first - 1];

    auto solver = [&](long long lambda){
        pair<long long, int> best = {LLONG_MAX, 0};
        for (int c = first; c <= n; c++) best = min(best, {f[c] + lambda * c, c});
        return best;
    };

    long long lo = 0, hi = 0;
    if (!slopes.empty()) lo = -slopes.back(), hi = -slopes.front();
    if (stress::rand_int(0, 1)) lo -= stress::rand_int(0, 1000000000000LL), hi += stress::rand_int(0, 1000000000000LL);

    for (int k = first; k <= n; k++) assert(aliens_trick(k, lo, hi, solver) == f[k]);
}

/// Split into k segments minimizing the sum of squared segment sums, against the O(n^2 k) DP
void check_squared_split(int n, int value_max){
    vector<long long> a(n), pre(n + 1, 0);
    for (auto& x : a) x = stress::rand_int(0, value_max);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + a[i];

    vector<vector<long long>> exact(n + 1, vector<long long>(n + 1, INF));
    exact[0][0] = 0;
    for (int c = 1; c <= n; c++){
        for (int j = 1; j <= n; j++){
            for (int i = 0; i < j; i++){
                if (exact[c - 1][i] == INF) continue;
                long long s = pre[j] - pre[i];
                exact[c][j] = min(exact[c][j], exact[c - 1][i] + s * s);
            }
        }
    }

    auto solver = [&](long long lambda){
        vector<pair<long long, int>> dp(n + 1, {LLONG_MAX, 0});
        dp[0] = {0, 0};
        for (int j = 1; j <= n; j++){
            for (int i = 0; i < j; i++){
                long long s = pre[j] - pre[i];
                dp[j] = min(dp[j], {dp[i].first + s * s + lambda, dp[i].second + 1});
            }
        }
        return dp[n];
    };

    for (int k = 1; k <= n; k++) assert(aliens_trick(k, 0LL, pre[n] * pre[n], solver) == exact[k][n]);
}

/// Maximum total of exactly k non-overlapping non-empty subarrays, against the O(n^2 k) DP
void check_max_k_subarrays(int n, int value_max){
    vector<long long> a(n), pre(n + 1, 0);
    for (auto& x : a) x = stress::rand_int(-value_max, value_max);
    for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + a[i];

    vector<vector<long long>> exact(n + 1, vector<long long>(n + 1, -INF));
    for (int j = 0; j <= n; j++) exact[0][j] = 0;
    for (int c = 1; c <= n; c++){
        for (int j = 1; j <= n; j++){
            exact[c][j] = exact[c][j - 1];
            for (int i = 0; i < j; i++){
                if (exact[c - 1][i] != -INF) exact[c][j] = max(exact[c][j], exact[c - 1][i] + pre[j] - pre[i]);
            }
        }
    }

    long long bound = 1;
    for (long long x : a) bound += abs(x);

    auto solver = [&](long long lambda){
        pair<long long, int> out = {0, 0}, in = {INF, 0};
        for (long long x : a){
            in = min(in, {out.first + lambda, out.second + 1});
            in.first -= x;
            out = min(out, in);
        }
        return out;
    };

    for (int k = 1; k <= n; k++) assert(-aliens_trick(k, -bound, bound, solver) == exact[k][n]);
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++) check_synthetic();

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 12), value_max = it % 3 == 0 ? 2 : (it % 3 == 1 ? 10 : 1000);
        check_squared_split(n, value_max);
        check_max_k_subarrays(n, value_max);
    }

    for (int it = 0; it < 3; it++){
        check_squared_split(70, 1000000);
        check_max_k_subarrays(70, 1000000000);
    }

    return 0;
}
