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
auto mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Get get){
    int q = queries.size();
    vector<decay_t<decltype(get(0))>> res(q);
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
        auto [l, r] = queries[qi];
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
    /// Distinct values per range: the callbacks keep a frequency table and a running count
    vector<int> a = {1, 2, 1, 3, 2, 2, 4};
    vector<pair<int, int>> queries = {{0, 6}, {0, 2}, {3, 5}, {4, 3}};
    vector<int> freq(5, 0);
    int count = 0;
    auto distinct = mo((int)a.size(), queries,
        [&](int i){ count += freq[a[i]]++ == 0; },
        [&](int i){ count -= --freq[a[i]] == 0; },
        [&](int){ return count; });
    assert((distinct == vector<int>{4, 2, 2, 0}));  /// {4, 3} is empty
    return 0;
}
