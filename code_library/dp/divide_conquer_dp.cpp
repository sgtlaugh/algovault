/***
 *
 * Divide and Conquer DP Optimization
 * Splits the prefix [0, i) into exactly k non-empty consecutive segments of minimum total cost
 *
 * Complexity: O(k n log n) calls to cost, O(n) memory
 *
 * divide_conquer_dp<T>(n, k, cost): dp[0..n] with dp[i] the cheapest split of the first i elements, k >= 0
 *
 * Notes:
 *   - cost(j, i) is the cost of the segment of elements j..i-1, called only for 0 <= j < i <= n
 *   - dp[i] = numeric_limits<T>::max() when i < k (fewer elements than segments)
 *   - Requires cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c) for a <= b <= c <= d, so the optimal j grows with i
 *   - Such costs: (sum of a non-negative segment)^2, sum of non-negative pair weights inside it, convex in its length
 *   - Terms f(j) + g(i) of any sign may be added to any such cost
 *   - T is passed explicitly and fits the sum of k costs, like long long
 *   - For maximization, negate the cost and the answer
 *   - knuth_partition_dp has the same call and is faster when k log2(n) > n, about 3x at k = n
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Cost>
vector<T> divide_conquer_dp(int n, int k, const Cost& cost){
    assert(k >= 0);

    const T inf = numeric_limits<T>::max();
    vector<T> prev(n + 1, inf), cur;
    if (k > n) return prev;
    if (k == 0){
        prev[0] = 0;
        return prev;
    }

    /// fills cur[l..r] knowing their optimal j lie in [opt_l, opt_r], invariant opt_l < l keeps the scan non-empty
    auto solve = [&](auto&& self, int l, int r, int opt_l, int opt_r) -> void{
        if (l > r) return;

        int mid = (l + r) / 2, best = -1;
        T best_val = inf;
        for (int j = opt_l; j <= min(mid - 1, opt_r); j++){
            T val = prev[j] + cost(j, mid);
            if (best == -1 || val < best_val) best = j, best_val = val;
        }
        cur[mid] = best_val;

        self(self, l, mid - 1, opt_l, best);
        self(self, mid + 1, r, best, opt_r);
    };

    for (int i = 1; i <= n; i++) prev[i] = cost(0, i);
    for (int g = 2; g <= k; g++){
        /// layer g - 1 is finite exactly from index g - 1, so starting j there keeps inf out of every addition
        cur.assign(n + 1, inf);
        solve(solve, g, n, g - 1, n - 1);
        prev.swap(cur);
    }

    return prev;
}

int main(){
    const long long inf = numeric_limits<long long>::max();

    /// Split {1, 2, 3, 4, 5} minimizing the sum of squared segment sums
    vector<long long> pre = {0, 1, 3, 6, 10, 15};
    auto squares = [&](int j, int i){ return (pre[i] - pre[j]) * (pre[i] - pre[j]); };

    vector<long long> two = divide_conquer_dp<long long>(5, 2, squares);
    assert(two[5] == 117);  /// {1, 2, 3} + {4, 5} = 36 + 81
    assert(two[3] == 18);   /// prefix {1, 2, 3} as {1, 2} + {3} = 9 + 9
    assert(two[1] == inf);  /// one element cannot make two segments

    assert(divide_conquer_dp<long long>(5, 3, squares)[5] == 77);  /// {1, 2, 3} + {4} + {5}
    return 0;
}
