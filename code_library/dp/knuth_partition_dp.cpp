/***
 *
 * Knuth Partition DP
 * Splits the prefix [0, i) into exactly k non-empty consecutive segments of minimum total cost
 *
 * Complexity: O(n^2 + n k) calls to cost whatever k is, O(n) memory
 *
 * knuth_partition_dp<T>(n, k, cost): dp[0..n] with dp[i] the cheapest split of the first i elements, k >= 0
 *
 * Notes:
 *   - cost(j, i) is the cost of the segment of elements j..i-1, called only for 0 <= j < i <= n
 *   - dp[i] = numeric_limits<T>::max() when i < k (fewer elements than segments)
 *   - Requires cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c) for a <= b <= c <= d
 *   - Such costs: (sum of a non-negative segment)^2, sum of non-negative pair weights inside it, convex in its length
 *   - Terms f(j) + g(i) of any sign may be added to any such cost
 *   - T is passed explicitly and fits the sum of k costs, like long long
 *   - For maximization, negate the cost and the answer
 *   - Same call, contract and result as divide_conquer_dp, so either can replace the other
 *   - Faster than divide_conquer_dp when k log2(n) > n, about 3x at k = n
 *   - Far slower for smaller k, since even k = 2 costs O(n^2) calls
 *   - A cost moving a window [j, i) by pointers makes about 2 moves per call, as in divide_conquer_dp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Cost>
vector<T> knuth_partition_dp(int n, int k, const Cost& cost){
    assert(k >= 0);

    const T inf = numeric_limits<T>::max();
    vector<T> dp(n + 1, inf);
    if (k > n) return dp;
    if (k == 0){
        dp[0] = 0;
        return dp;
    }

    for (int i = n; i >= 1; i--) dp[i] = cost(0, i);

    /// opt[n + 1] = n - 1 caps the first scan of every layer, layer 1 has opt[i] = 0
    vector<int> opt(n + 2, 0);
    opt[n + 1] = n - 1;
    for (int g = 2; g <= k; g++){
        /// i falls, so dp[j < i] and opt[i] still hold layer g - 1 while opt[i + 1] already holds layer g
        for (int i = n; i >= g; i--){
            int best = max(opt[i], g - 1), hi = min(opt[i + 1], i - 1);
            T best_val = dp[best] + cost(best, i);
            for (int j = best + 1; j <= hi; j++){
                T val = dp[j] + cost(j, i);
                if (val < best_val) best = j, best_val = val;
            }

            dp[i] = best_val;
            opt[i] = best;
        }
        dp[g - 1] = inf;
    }

    return dp;
}

int main(){
    const long long inf = numeric_limits<long long>::max();

    /// Split {1, 2, 3, 4, 5} minimizing the sum of squared segment sums
    vector<long long> pre = {0, 1, 3, 6, 10, 15};
    auto squares = [&](int j, int i){ return (pre[i] - pre[j]) * (pre[i] - pre[j]); };

    vector<long long> two = knuth_partition_dp<long long>(5, 2, squares);
    assert(two[5] == 117);  /// {1, 2, 3} + {4, 5} = 36 + 81
    assert(two[3] == 18);   /// prefix {1, 2, 3} as {1, 2} + {3} = 9 + 9
    assert(two[1] == inf);  /// one element cannot make two segments

    assert(knuth_partition_dp<long long>(5, 3, squares)[5] == 77);  /// {1, 2, 3} + {4} + {5}
    assert(knuth_partition_dp<long long>(5, 5, squares)[5] == 55);  /// every element alone, 1 + 4 + 9 + 16 + 25
    return 0;
}
