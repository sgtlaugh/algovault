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
    int n, m;
    scanf("%d %d", &n, &m);

    vector<long long> pos(n, 0);
    for (int i = 1; i < n; i++){
        scanf("%lld", &pos[i]);
        pos[i] += pos[i - 1];
    }

    int lg = 1;
    while ((1 << lg) <= n) lg++;
    vector<vector<int>> st(lg, vector<int>(n * m));
    for (int i = 0; i < n * m; i++) scanf("%d", &st[0][i]);
    for (int h = 1; h < lg; h++){
        for (int i = 0; i + (1 << h) <= n; i++){
            for (int t = 0; t < m; t++){
                st[h][i * m + t] = max(st[h - 1][i * m + t], st[h - 1][(i + (1 << (h - 1))) * m + t]);
            }
        }
    }
    vector<int> log_of(n + 1, 0);
    for (int i = 2; i <= n; i++) log_of[i] = log_of[i / 2] + 1;

    /// a dummy element 0 makes layer 1 free, so layer 2 picks the best left end l = j - 1 for every right end r = i - 2
    /// interval max sums are inverse Monge, so the negated happiness satisfies the quadrangle inequality
    auto cost = [&](int j, int i) -> long long{
        if (j == 0) return 0;
        int l = j - 1, r = i - 2, h = log_of[r - l + 1];
        const int* a = &st[h][l * m];
        const int* b = &st[h][(r - (1 << h) + 1) * m];
        long long res = pos[r] - pos[l];
        for (int t = 0; t < m; t++) res -= max(a[t], b[t]);
        return res;
    };

    vector<long long> dp = divide_conquer_dp<long long>(n + 1, 2, cost);
    printf("%lld\n", -*min_element(dp.begin() + 2, dp.end()));
    return 0;
}
