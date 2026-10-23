/***
 *
 * Aliens Trick (WQS Binary Search, Lagrangian Relaxation)
 * Minimum cost of a solution with exactly k items, when that cost f(k) is convex in k
 *
 * Complexity: O(log(hi - lo)) solver calls
 *
 * aliens_trick(k, lo, hi, solver): f(k), in the solver's value type
 * solver(lambda): {min over all solutions of cost + lambda * items, the fewest items among solutions reaching it}
 *
 * Notes:
 *   - Integer costs, f convex on its domain and k inside it
 *   - [lo, hi] contains every -(f(c) - f(c - 1)), so M = max |f(c) - f(c - 1)| gives [-M, M]
 *   - The bounds' type fits hi - lo, the solver's value type fits cost + max(|lo|, |hi|) * (max items)
 *   - The count must be the fewest: several counts can tie for one lambda, so count == k may never come back
 *   - A DP over pairs {cost, items} compared lexicographically returns exactly that pair
 *   - Maximizing a concave f: minimize the negated gains and negate the result
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
