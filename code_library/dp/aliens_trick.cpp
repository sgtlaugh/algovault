/***
 *
 * Aliens Trick (WQS Binary Search, Lagrangian Relaxation)
 * Minimum cost of a solution with exactly k items, when that cost f(k) is convex in k
 *
 * Complexity: O(log(hi - lo)) calls to the solver, plus one
 *
 * aliens_trick(k, lo, hi, solver) returns f(k)
 *   solver(lambda) returns {value, count}: value = min over all solutions of cost + lambda * items,
 *   count = the MINIMUM number of items among the solutions reaching that value
 *   (a DP over pairs {cost, items} compared lexicographically gives exactly this)
 *
 * Requirements:
 *   - integer costs, f convex on its domain and k inside the domain
 *   - [lo, hi] contains every marginal -(f(c) - f(c - 1)), so M = max |f(c) - f(c - 1)| gives [-M, M]
 *   - the bounds' type is wide enough for hi - lo, the solver's value type for cost + max(|lo|, |hi|) * (max items)
 *   - f(k) comes back in the solver's value type, so int bounds with an int64 solver still return int64
 *
 * Ties: several counts can be optimal for one lambda, so a lambda with count == k may not exist
 * The driver takes the smallest lambda whose minimum optimal count is <= k; k is then optimal there too,
 * and f(k) = value - lambda * k (never value - lambda * count)
 *
 * Maximizing a concave f: minimize the negated gains and negate the result
 *
 * Example: split a non-negative array into k segments minimizing the sum of squared segment sums
 *   aliens_trick(k, 0LL, S * S, solver), solver(lambda) runs the O(n^2) split DP with + lambda per segment
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Solver>
auto aliens_trick(int k, T lo, T hi, Solver solver){
    while (lo < hi){
        T mid = lo + (hi - lo) / 2;
        if (solver(mid).second <= k) hi = mid;
        else lo = mid + 1;
    }

    auto value = solver(lo).first;
    return value - decltype(value)(lo) * k;
}

int main(){
    auto min_squared_split = [](const vector<long long>& a, int k){
        int n = a.size();
        vector<long long> pre(n + 1, 0);
        for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + a[i];

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

        return aliens_trick(k, 0LL, pre[n] * pre[n], solver);
    };

    auto max_k_subarrays = [](const vector<long long>& a, int k){
        long long bound = 1;
        for (long long x : a) bound += abs(x);

        auto solver = [&](long long lambda){
            pair<long long, int> out = {0, 0}, in = {LLONG_MAX / 2, 0};
            for (long long x : a){
                in = min(in, {out.first + lambda, out.second + 1});
                in.first -= x;
                out = min(out, in);
            }
            return out;
        };

        return -aliens_trick(k, -bound, bound, solver);
    };

    assert(min_squared_split({1, 2, 3, 4}, 2) == 52);  /// (1 + 2 + 3)^2 + 4^2
    assert(min_squared_split({1, 2, 3, 4}, 3) == 34);  /// 3^2 + 3^2 + 4^2
    assert(min_squared_split({2, 0, 0, 2}, 3) == 8);   /// f(2) = f(3), no lambda has minimum count 3

    assert(max_k_subarrays({3, -1, 2, -5, 4}, 2) == 8);  /// [3, -1, 2] and [4]
    assert(max_k_subarrays({3, -1, 2, -5, 4}, 3) == 9);  /// [3], [2] and [4]
    assert(max_k_subarrays({-3, -1, -2}, 2) == -3);      /// exactly k, even when every pick loses

    auto shifted_squares = [](long long lambda){      /// f(c) = 5e9 + c^2 for c in [0, 10]
        pair<long long, int> best = {LLONG_MAX, 0};
        for (int c = 0; c <= 10; c++) best = min(best, {5000000000LL + c * c + lambda * c, c});
        return best;
    };
    assert(aliens_trick(2, -20, 20, shifted_squares) == 5000000004LL);  /// int bounds keep the int64 value

    auto steep = [](long long lambda){                /// f(c) = c^2 - 1.2e9 c for c in [0, 10]
        pair<long long, int> best = {LLONG_MAX, 0};
        for (int c = 0; c <= 10; c++) best = min(best, {1LL * c * c - 1200000000LL * c + lambda * c, c});
        return best;
    };
    assert(aliens_trick(5, 1000000000, 1300000000, steep) == -5999999975LL);  /// lambda * k overflows int

    return 0;
}
