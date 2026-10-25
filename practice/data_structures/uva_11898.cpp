/***
 *
 * Mo's Algorithm
 * Answers offline range queries on an array by moving a window over it
 * Variants: mo_with_updates.cpp (point assignments), mo_with_rollback.cpp (add-only state), tree_path_mo.cpp (tree paths)
 *
 * Complexity: O(n * sqrt(q) + q) calls to add / remove, plus O(q log q) sorting
 *
 * mo(n, queries, add, remove, get) returns the answers, res[qi] = get(qi) taken while the window is exactly query qi:
 *   queries[i] = {l, r}, 0-based inclusive ranges over positions 0..n-1, r = l - 1 is an empty range
 *   add(i) / remove(i) move position i into / out of the window, get(qi) returns the answer for the current window
 *   The window starts and ends empty, so the same callbacks and state can serve several runs
 *
 *   auto distinct = mo(n, queries, [&](int i){ count += freq[a[i]]++ == 0; },
 *                      [&](int i){ count -= --freq[a[i]] == 0; }, [&](int){ return count; });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED mo
template <typename Add, typename Remove, typename Get>
auto mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Get get) -> vector<typename decay<decltype(get(0))>::type>{
    int q = queries.size();
    vector<typename decay<decltype(get(0))>::type> res(q);
    if (q == 0) return res;
    int block = max(1, (int)(n / sqrt((double)q)));

    /// Odd-even block order: Hilbert curve order measured about 5% more add / remove calls on random n = q = 5e5
    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? queries[a].second > queries[b].second : queries[a].second < queries[b].second;
    });

    /// Growing before shrinking keeps every removed position inside the window, empty targets included
    int cur_l = 0, cur_r = -1;
    for (int qi : order){
        int l = queries[qi].first, r = queries[qi].second;
        assert(0 <= l && l <= r + 1 && r < n);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        res[qi] = get(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);  /// leave the window empty so the caller's state can be reused
    return res;
}
/// END SHARED mo

int main(){
    const int MAXV = 10000;
    int t;

    scanf("%d", &t);
    while (t--){
        int n, q;
        scanf("%d", &n);
        vector<int> a(n);
        for (auto& x : a) scanf("%d", &x);

        scanf("%d", &q);
        vector<pair<int, int>> queries(q);
        for (auto& e : queries){
            scanf("%d %d", &e.first, &e.second);
            e.first--, e.second--;
        }

        vector<int> cnt(MAXV + 1, 0);
        int repeated = 0;
        auto res = mo(n, queries, [&](int i){ repeated += ++cnt[a[i]] == 2; },
                      [&](int i){ repeated -= cnt[a[i]]-- == 2; },
                      [&](int){
                          if (repeated) return 0;
                          int best = MAXV, last = 0;
                          for (int v = 1; v <= MAXV; v++){
                              if (!cnt[v]) continue;
                              if (last) best = min(best, v - last);
                              last = v;
                          }
                          return best;
                      });
        for (int x : res) printf("%d\n", x);
    }
    return 0;
}
