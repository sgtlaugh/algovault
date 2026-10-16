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
    int n, m, k;
    scanf("%d %d %d", &n, &m, &k);
    vector<int> a(n);
    for (int& v : a) scanf("%d", &v);
    vector<pair<int, int>> queries(k);
    for (auto& [l, r] : queries){
        scanf("%d %d", &l, &r);
        l--, r--;
    }

    /// The window only grows, so a value's span in it is its first and last position, with first = -1 when absent
    struct Change{
        int value, first, last, best;
    };
    vector<int> first(m + 1, -1), last(m + 1, -1);
    vector<Change> history;
    vector<size_t> marks;
    int best = 0;
    auto res = mo_with_rollback(n, queries,
        [&](int i){
            int v = a[i];
            history.push_back({v, first[v], last[v], best});
            if (first[v] == -1) first[v] = last[v] = i;
            else first[v] = min(first[v], i), last[v] = max(last[v], i);
            best = max(best, last[v] - first[v]);
        },
        [&](){ marks.push_back(history.size()); },
        [&](){
            while (history.size() > marks.back()){
                Change c = history.back();
                history.pop_back();
                first[c.value] = c.first, last[c.value] = c.last, best = c.best;
            }
            marks.pop_back();
        },
        [&](int){ return best; });
    for (int v : res) printf("%d\n", v);
    return 0;
}
