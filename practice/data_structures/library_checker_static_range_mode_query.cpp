/***
 *
 * Mo's Algorithm With Rollback
 * Answers offline range queries when the state can add a position but not remove one, like a max or a DSU with rollback
 *
 * Complexity: O(n sqrt(q) + q) calls to add plus O(q log q) sorting, a rollback costs as much as undoing its adds
 *
 * mo_with_rollback(n, queries, add, snapshot, rollback, get): the answers, res[qi] = get(qi) over the window queries[qi]
 * queries[qi] = {l, r}: the inclusive range l..r of positions 0..n - 1, r = l - 1 is an empty range
 * add(i): puts position i into the window, which stays contiguous and grows at one end
 * snapshot(): saves the state on a stack, at most 2 are saved at once
 * rollback(): restores the top saved state and pops it
 * get(qi): the answer for the current window
 *
 * The window starts and ends empty, and empty ranges are answered first, while the state is still empty
 * A history stack gives snapshot and rollback for any add, main() shows it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Add, typename Snapshot, typename Rollback, typename Get>
auto mo_with_rollback(int n, const vector<pair<int, int>>& queries, Add add, Snapshot snapshot, Rollback rollback, Get get){
    int q = queries.size();
    vector<decay_t<decltype(get(0))>> res(q);
    vector<int> order;
    for (int qi = 0; qi < q; qi++){
        auto [l, r] = queries[qi];
        assert(0 <= l && l <= r + 1 && r < n);
        if (l > r) res[qi] = get(qi);
        else order.push_back(qi);
    }

    int m = order.size();
    if (m == 0) return res;
    int block = max(1, (int)(n / sqrt((double)m)));
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return queries[a].second < queries[b].second;
    });

    /// Per block of l the right part grows from the block end and is kept, the left part is added and rolled back per query
    for (int i = 0; i < m; ){
        int b = queries[order[i]].first / block, block_end = min(n - 1, (b + 1) * block - 1), cur_r = block_end;
        snapshot();
        for (; i < m && queries[order[i]].first / block == b; i++){
            int qi = order[i];
            auto [l, r] = queries[qi];
            while (cur_r < r) add(++cur_r);
            snapshot();
            if (r <= block_end){
                for (int j = l; j <= r; j++) add(j);  /// sorted by r, so these run before the right part grows
            }
            else{
                for (int j = block_end; j >= l; j--) add(j);
            }
            res[qi] = get(qi);
            rollback();
        }
        rollback();
    }

    return res;
}

int main(){
    int n, q;
    scanf("%d %d", &n, &q);
    vector<int> a(n);
    for (int& v : a) scanf("%d", &v);
    vector<int> sorted_a = a;
    sort(sorted_a.begin(), sorted_a.end());
    for (int& v : a) v = lower_bound(sorted_a.begin(), sorted_a.end(), v) - sorted_a.begin();

    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries){
        scanf("%d %d", &l, &r);
        r--;
    }

    /// A mode can only change on an add, to the added value, so the history keeps the old best per add
    vector<int> freq(n, 0);
    vector<pair<int, pair<int, int>>> history;
    vector<size_t> marks;
    pair<int, int> best = {0, 0};
    auto res = mo_with_rollback(n, queries,
        [&](int i){
            history.push_back({a[i], best});
            if (++freq[a[i]] > best.second) best = {a[i], freq[a[i]]};
        },
        [&](){ marks.push_back(history.size()); },
        [&](){
            while (history.size() > marks.back()){
                freq[history.back().first]--;
                best = history.back().second;
                history.pop_back();
            }
            marks.pop_back();
        },
        [&](int){ return best; });
    for (auto [mode, count] : res) printf("%d %d\n", sorted_a[mode], count);
    return 0;
}
