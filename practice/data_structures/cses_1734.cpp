/***
 *
 * Mo's Algorithm
 * Answers offline range queries on an array by moving a window over it
 *
 * Complexity: O(n sqrt(q) + q) calls to add and remove, plus O(q log q) sorting
 *
 * mo(n, queries, add, remove, get): the answers, res[i] = get(i) taken while the window holds exactly queries[i] = {l, r}, 0 <= l <= r + 1 <= n
 *
 * add(i) and remove(i) move position i into and out of the window, get(i) returns the answer for the current window
 * r = l - 1 is an empty range
 * The window starts and ends empty, so the same callbacks and state can serve several runs
 * Use mo_with_updates.cpp for point assignments, mo_with_rollback.cpp for add-only state, tree_path_mo.cpp for tree paths
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
    int n, q;
    scanf("%d %d", &n, &q);
    vector<int> x(n);
    for (int& v : x) scanf("%d", &v);
    vector<int> sorted_x = x;
    sort(sorted_x.begin(), sorted_x.end());
    for (int& v : x) v = lower_bound(sorted_x.begin(), sorted_x.end(), v) - sorted_x.begin();

    vector<pair<int, int>> queries(q);
    for (auto& [a, b] : queries){
        scanf("%d %d", &a, &b);
        a--, b--;
    }

    vector<int> freq(n, 0);
    int distinct = 0;
    auto res = mo(n, queries,
        [&](int i){ distinct += freq[x[i]]++ == 0; },
        [&](int i){ distinct -= --freq[x[i]] == 0; },
        [&](int){ return distinct; });
    for (int v : res) printf("%d\n", v);
    return 0;
}
