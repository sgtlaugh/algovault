/***
 *
 * Knuth Optimization
 * Interval DP dp[i][j] = min over i < k < j of dp[i][k] + dp[k][j] + cost(i, j), dp[i][i + 1] = 0
 *
 * Complexity: O(n^2) time and memory
 *
 * knuth_dp<T>(n, cost): dp table over the points 0..n, dp[i][j] for 0 <= i < j <= n, the answer is dp[0][n]
 * cost(i, j) is called once per interval with 0 <= i, i + 2 <= j <= n and returns T
 *
 * Valid when for all a <= b <= c <= d:
 *     cost(b, c) <= cost(a, d)                                  (monotone on inclusion)
 *     cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c)        (quadrangle inequality)
 * Then the smallest optimal split opt[i][j] satisfies opt[i][j - 1] <= opt[i][j] <= opt[i + 1][j],
 * so each diagonal of the table scans O(n) splits in total
 *
 * Any f(sum of a[i..j - 1]) with f convex, nondecreasing and a >= 0 qualifies, interval sums included
 * A leaf cost on every [i, i + 1] adds the sum of the leaf costs over i..j - 1 to dp[i][j] regardless of the split tree,
 * so add it after the DP
 *
 * Merging n adjacent piles with cost = size of the merged pile:
 *     auto dp = knuth_dp<long long>(n, [&](int i, int j){ return pre[j] - pre[i]; });
 * Optimal BST over keys 1..m with frequencies f: points 0..m + 1, cost(i, j) = f[i + 1] + ... + f[j - 1]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Cost>
vector<vector<T>> knuth_dp(int n, Cost cost){
    vector<vector<T>> dp(n + 1, vector<T>(n + 1, 0));
    vector<vector<int>> opt(n + 1, vector<int>(n + 1, 0));
    for (int i = 0; i < n; i++) opt[i][i + 1] = i;

    for (int len = 2; len <= n; len++){
        for (int i = 0, j = len; j <= n; i++, j++){
            int best = max(i + 1, opt[i][j - 1]);
            T best_val = dp[i][best] + dp[best][j];
            for (int k = best + 1; k <= opt[i + 1][j]; k++){
                T val = dp[i][k] + dp[k][j];
                if (val < best_val) best_val = val, best = k;
            }

            dp[i][j] = best_val + cost(i, j);
            opt[i][j] = best;
        }
    }

    return dp;
}

int main(){
    auto merge_piles = [](const vector<long long>& piles){
        vector<long long> pre(piles.size() + 1, 0);
        for (size_t i = 0; i < piles.size(); i++) pre[i + 1] = pre[i] + piles[i];
        return knuth_dp<long long>(piles.size(), [&](int i, int j){ return pre[j] - pre[i]; });
    };

    auto dp = merge_piles({1, 2, 3});
    assert(dp[0][3] == 9 && dp[0][2] == 3 && dp[1][3] == 5 && dp[0][1] == 0);
    assert(merge_piles({4, 1, 1, 4})[0][4] == 18);
    assert(merge_piles({7})[0][1] == 0);
    assert(merge_piles({})[0][0] == 0);
    assert(merge_piles(vector<long long>(2048, 1))[0][2048] == 2048 * 11);

    vector<long long> freq_pre = {0, 34, 42, 92};
    assert((knuth_dp<long long>(4, [&](int i, int j){ return freq_pre[j - 1] - freq_pre[i]; })[0][4] == 142));

    return 0;
}
