/***
 *
 * Mo's Algorithm With Updates
 * Answers offline range queries on an array that also receives point assignments, by moving a window in time as well
 *
 * Complexity: O((q + u) * n^(2/3) + n^(4/3)) calls to add / remove for q queries and u updates, plus O(q log q) sorting
 *
 * mo_with_updates(a, queries, updates, add, remove, get) returns the answers, res[qi] = get(qi):
 *   queries[i] = {l, r, t}: the inclusive range l..r of a after the first t updates, 0 <= t <= u, r = l - 1 is empty
 *   updates[j] = {pos, value}: update j assigns a[pos] = value
 *   add(i) / remove(i) move position i into / out of the window and read the current a[i]
 *   get(qi) returns the answer for the current window
 *   Updates are applied and undone by swapping values in a, inside the window wrapped in remove / add
 *   The run starts and ends with the window empty and no update applied, so a is left as it was
 *   The n^(2/3) block assumes q and u close to n, with no updates plain mo.cpp is faster: n * sqrt(q) against q * n^(2/3)
 *
 *   auto freq_of_x = mo_with_updates(a, queries, updates, [&](int i){ cnt[a[i]]++; },
 *                                    [&](int i){ cnt[a[i]]--; }, [&](int qi){ return cnt[x[qi]]; });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Add, typename Remove, typename Get>
auto mo_with_updates(vector<T>& a, const vector<array<int, 3>>& queries, const vector<pair<int, T>>& updates, Add add, Remove remove, Get get){
    int n = a.size(), q = queries.size(), u = updates.size();
    vector<decay_t<decltype(get(0))>> res(q);
    if (q == 0) return res;
    int block = max(1, (int)cbrt((double)n * n));

    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y){
        int lx = queries[x][0] / block, ly = queries[y][0] / block;
        if (lx != ly) return lx < ly;
        int rx = queries[x][1] / block, ry = queries[y][1] / block;
        if (rx != ry) return rx < ry;
        return queries[x][2] < queries[y][2];
    });

    /// Updates are applied in increasing and undone in decreasing order, so swapping with the stored value does both
    vector<T> value(u);
    for (int j = 0; j < u; j++) value[j] = updates[j].second;

    int cur_l = 0, cur_r = -1, cur_t = 0;
    auto flip = [&](int j){
        int p = updates[j].first;
        assert(0 <= p && p < n);
        bool inside = cur_l <= p && p <= cur_r;
        if (inside) remove(p);
        swap(a[p], value[j]);
        if (inside) add(p);
    };

    for (int qi : order){
        auto [l, r, t] = queries[qi];
        assert(0 <= l && l <= r + 1 && r < n && 0 <= t && t <= u);
        while (cur_t < t) flip(cur_t++);
        while (cur_t > t) flip(--cur_t);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        res[qi] = get(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);
    while (cur_t > 0) flip(--cur_t);
    return res;
}

int main(){
    /// Distinct values per range: after update 0, b = {1, 1, 1, 3, 2}
    vector<int> b = {1, 2, 1, 3, 2}, freq(5, 0);
    vector<pair<int, int>> updates = {{1, 1}};
    vector<array<int, 3>> queries = {{0, 2, 0}, {0, 2, 1}, {3, 2, 1}};
    int count = 0;
    auto distinct = mo_with_updates(b, queries, updates,
        [&](int i){ count += freq[b[i]]++ == 0; },
        [&](int i){ count -= --freq[b[i]] == 0; },
        [&](int){ return count; });
    assert((distinct == vector<int>{2, 1, 0}));
    assert((b == vector<int>{1, 2, 1, 3, 2}));  /// every update is undone at the end

    return 0;
}
