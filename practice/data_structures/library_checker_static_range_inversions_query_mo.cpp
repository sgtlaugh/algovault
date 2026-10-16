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
    vector<int> a(n);
    for (int& v : a) scanf("%d", &v);
    vector<int> sorted_a = a;
    sort(sorted_a.begin(), sorted_a.end());
    for (int& v : a) v = lower_bound(sorted_a.begin(), sorted_a.end(), v) - sorted_a.begin() + 1;

    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        scanf("%d %d", &l, &r);
        r--;
    }

    vector<int> fenwick(n + 1, 0);
    auto update = [&](int p, int d){
        for (; p <= n; p += p & -p) fenwick[p] += d;
    };
    auto prefix = [&](int p){
        int s = 0;
        for (; p > 0; p -= p & -p) s += fenwick[p];
        return s;
    };

    /// add and remove only get a position, so the window bounds are tracked here to know which end moves:
    /// an element entering on the left inverts with the smaller ones, on the right with the larger ones
    int lo = 0, hi = -1;
    long long inversions = 0;
    auto res = mo(n, queries,
        [&](int i){
            if (lo > hi) lo = hi = i;
            else if (i < lo) lo = i, inversions += prefix(a[i] - 1);
            else hi = i, inversions += (hi - lo) - prefix(a[i]);
            update(a[i], 1);
        },
        [&](int i){
            update(a[i], -1);
            if (lo == hi) lo = 0, hi = -1;
            else if (i == lo) lo++, inversions -= prefix(a[i] - 1);
            else hi--, inversions -= (hi - lo + 1) - prefix(a[i]);
        },
        [&](int){ return inversions; });
    for (long long v : res) printf("%lld\n", v);
    return 0;
}
