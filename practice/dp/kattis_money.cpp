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
    int m, n;
    scanf("%d %d", &m, &n);
    vector<pair<long long, long long>> prod(m), cons(n);
    for (auto& x: prod) scanf("%lld %lld", &x.second, &x.first);
    for (auto& x: cons) scanf("%lld %lld", &x.second, &x.first);

    /// keep only undominated companies: producers by day ascending with falling price, consumers by end ascending with falling price
    sort(prod.begin(), prod.end());
    vector<long long> pd, pp;
    for (auto& x: prod){
        if (pp.empty() || x.second < pp.back()) pd.push_back(x.first), pp.push_back(x.second);
    }
    sort(cons.rbegin(), cons.rend());
    vector<long long> ce, cq;
    for (auto& x: cons){
        if (!cq.empty() && x.second <= cq.back()) continue;
        /// a consumer without any producer starting earlier and selling cheaper would make the raw product lie
        int idx = lower_bound(pd.begin(), pd.end(), x.first) - pd.begin() - 1;
        if (idx >= 0 && pp[idx] < x.second) ce.push_back(x.first), cq.push_back(x.second);
    }
    reverse(ce.begin(), ce.end());
    reverse(cq.begin(), cq.end());
    if (ce.empty()){
        puts("0");
        return 0;
    }

    /// indices 1..P are producers and P+1..P+C consumers; every index gets both roles with monotone sentinel attributes,
    /// so (q - p)(e - d) stays inverse Monge on the whole range and a consumer index used as a producer is hopeless
    const long long far = 2000000000LL;
    int np = pd.size(), nc = ce.size();
    auto cost = [&](int j, int i) -> long long{
        if (j == 0) return 0;
        long long p = j <= np ? pp[j - 1] : pp.back(), d = j <= np ? pd[j - 1] : far;
        long long q = i <= np ? far : cq[i - np - 1], e = i <= np ? 0 : ce[i - np - 1];
        return -(q - p) * (e - d);
    };

    vector<long long> dp = divide_conquer_dp<long long>(np + nc, 2, cost);
    printf("%lld\n", max(0LL, -*min_element(dp.begin() + np + 1, dp.end())));
    return 0;
}
