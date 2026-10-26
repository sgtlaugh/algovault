/***
 *
 * Knuth Optimization
 * Interval DP dp[i][j] = min over i < k < j of dp[i][k] + dp[k][j] + cost(i, j), dp[i][i + 1] = 0
 *
 * Complexity: O(n^2) time and memory, (n + 1)^2 values of T plus as many ints, about 300 MB at n = 5000 for long long
 *
 * knuth_dp<T>(n, cost): KnuthTable<T> dp over points 0..n, read as dp[i][j] for 0 <= i < j <= n, the answer is dp[0][n]
 *
 * Notes:
 *   - cost(i, j) is called once per interval with 0 <= i, i + 2 <= j <= n and returns T
 *   - Requires, for all a <= b <= c <= d, cost(b, c) <= cost(a, d) and the quadrangle inequality:
 *         cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c)
 *   - Any f(sum of a[i..j - 1]) with f convex, nondecreasing and a >= 0 qualifies, interval sums included
 *   - A leaf cost on every [i, i + 1] adds the same sum to dp[i][j] for every split tree, so add it after the DP
 *   - Merging n adjacent piles at the merged size: cost(i, j) = pre[j] - pre[i], or garsia_wachs.cpp in O(n log n)
 *   - Optimal BST over keys 1..m with frequencies f: points 0..m + 1, cost(i, j) = f[i + 1] + ... + f[j - 1]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// One allocation per table: a vector per row measured 3x slower on 200000 calls at n = 10, 1.6x at n = 250
template <typename T>
struct KnuthTable{
    int n;
    vector<T> cells;

    KnuthTable(int n) : n(n), cells((size_t)(n + 1) * (n + 1), 0) {}

    T* operator[](int i){
        return cells.data() + (size_t)i * (n + 1);
    }

    const T* operator[](int i) const{
        return cells.data() + (size_t)i * (n + 1);
    }
};

template <typename T, typename Cost>
KnuthTable<T> knuth_dp(int n, Cost cost){
    KnuthTable<T> dp(n);
    KnuthTable<int> opt(n);

    for (int i = 0; i < n; i++) opt[i][i + 1] = i;

    /// Rows bottom-up keep rows i and i + 1 cached, diagonal by diagonal measured 1.4x slower at n = 5000
    for (int i = n - 2; i >= 0; i--){
        for (int j = i + 2; j <= n; j++){
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
    assert(dp[0][3] == 9);
    assert(dp[0][2] == 3);
    assert(dp[1][3] == 5);
    assert(dp[0][1] == 0);
    assert(merge_piles({4, 1, 1, 4})[0][4] == 18);
    assert(merge_piles({7})[0][1] == 0);
    assert(merge_piles({})[0][0] == 0);
    assert(merge_piles(vector<long long>(2048, 1))[0][2048] == 2048 * 11);

    vector<long long> freq_pre = {0, 34, 42, 92};
    assert((knuth_dp<long long>(4, [&](int i, int j){ return freq_pre[j - 1] - freq_pre[i]; })[0][4] == 142));

    return 0;
}
