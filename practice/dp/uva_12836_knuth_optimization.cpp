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

vector<int> lis_ending_at(const vector<int>& a){
    vector<int> tail, len;

    for (size_t i = 0; i < a.size(); i++){
        int pos = lower_bound(tail.begin(), tail.end(), a[i]) - tail.begin();
        if (pos == (int)tail.size()) tail.push_back(a[i]);
        else tail[pos] = a[i];
        len.push_back(pos + 1);
    }

    return len;
}

int main(){
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);
        vector<int> a(n);
        for (int i = 0; i < n; i++) scanf("%d", &a[i]);

        vector<int> left = lis_ending_at(a);
        vector<int> right = lis_ending_at(vector<int>(a.rbegin(), a.rend()));

        vector<long long> pre(n + 1, 0);
        for (int i = 0; i < n; i++) pre[i + 1] = pre[i] + left[i] + right[n - 1 - i] - 1;

        auto dp = knuth_dp<long long>(n, [&](int i, int j){ return pre[j] - pre[i]; });
        printf("Case %d: %lld\n", tc, dp[0][n]);
    }

    return 0;
}
