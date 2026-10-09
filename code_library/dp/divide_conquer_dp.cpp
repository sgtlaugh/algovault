/***
 *
 * Divide and Conquer DP Optimization
 * Splits the prefix [0, i) into exactly k non-empty consecutive segments of minimum total cost
 *
 * Complexity: O(k * n log n) calls to cost, O(n) memory
 *
 * divide_conquer_dp<T>(n, k, cost) returns dp of size n + 1 for 0 <= k, where
 *   dp[i] = min over j < i of (layer k - 1 at j) + cost(j, i), layer 0 being 0 at i = 0
 *   cost(j, i) is the cost of the segment of elements j..i-1, called only for 0 <= j < i <= n
 *   dp[i] = numeric_limits<T>::max() when i < k (fewer elements than segments)
 *
 * Correct only when the optimal j is non-decreasing in i, which holds when cost satisfies
 * the quadrangle inequality cost(a, c) + cost(b, d) <= cost(a, d) + cost(b, c) for a <= b <= c <= d
 * Common such costs: (sum of a non-negative segment)^2, sum of non-negative pair weights inside the segment,
 * any convex function of the segment length, plus terms f(j) + g(i) of any sign
 * Pass T explicitly and wide enough for the sum of k costs, e.g. long long
 * For maximization, negate the cost and the answer
 *
 * Example: split {1, 2, 3, 4, 5} into 2 segments minimizing the sum of squared segment sums
 *   vector<long long> pre = {0, 1, 3, 6, 10, 15};
 *   auto dp = divide_conquer_dp<long long>(5, 2, [&](int j, int i){ return (pre[i] - pre[j]) * (pre[i] - pre[j]); });
 *   // dp[5] = 117 from {1, 2, 3} + {4, 5}
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
    vector<long long> pre = {0, 1, 3, 6, 10, 15};
    auto squares = [&](int j, int i){ return (pre[i] - pre[j]) * (pre[i] - pre[j]); };

    assert((divide_conquer_dp<long long>(5, 1, squares) == vector<long long>{inf, 1, 9, 36, 100, 225}));
    assert((divide_conquer_dp<long long>(5, 2, squares) == vector<long long>{inf, inf, 5, 18, 52, 117}));
    assert(divide_conquer_dp<long long>(5, 3, squares)[5] == 77);
    assert(divide_conquer_dp<long long>(5, 5, squares)[5] == 55);
    assert((divide_conquer_dp<long long>(5, 0, squares) == vector<long long>{0, inf, inf, inf, inf, inf}));
    assert((divide_conquer_dp<long long>(5, 6, squares) == vector<long long>(6, inf)));
    assert((divide_conquer_dp<long long>(0, 0, squares) == vector<long long>{0}));
    assert((divide_conquer_dp<long long>(0, 1, squares) == vector<long long>{inf}));

    /// Codeforces 321E Ciel and Gondolas samples: cost is the unfamiliarity sum over pairs inside a gondola
    auto gondolas = [](const vector<vector<int>>& u, int k){
        int n = u.size();
        vector<vector<long long>> sum(n + 1, vector<long long>(n + 1, 0));
        for (int a = 0; a < n; a++){
            for (int b = 0; b < n; b++) sum[a + 1][b + 1] = sum[a][b + 1] + sum[a + 1][b] - sum[a][b] + u[a][b];
        }

        auto cost = [&](int j, int i){ return (sum[i][i] - sum[j][i] - sum[i][j] + sum[j][j]) / 2; };
        return divide_conquer_dp<long long>(n, k, cost)[n];
    };

    vector<vector<int>> strangers(8, vector<int>(8, 1));
    for (int a = 0; a < 8; a++) strangers[a][a] = 0;

    assert(gondolas({{0, 0, 1, 1, 1}, {0, 0, 1, 1, 1}, {1, 1, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 1, 0, 0, 0}}, 2) == 0);
    assert(gondolas(strangers, 3) == 7);
    assert(gondolas({{0, 2, 0}, {2, 0, 3}, {0, 3, 0}}, 2) == 2);

    return 0;
}
