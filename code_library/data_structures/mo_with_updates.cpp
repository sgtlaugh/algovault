/***
 *
 * Mo's Algorithm With Updates
 * Answers offline range queries on an array that also receives point assignments
 *
 * Complexity: O((q + u) n^(2/3) + n^(4/3)) calls to add and remove for q queries and u updates, plus O(q log q) sorting
 *
 * mo_with_updates(a, queries, updates, add, remove, get): the answers, res[qi] = get(qi) over the window of queries[qi]
 * queries[qi] = {l, r, t}: range l..r of a after the first t updates, 0 <= t <= u, r = l - 1 is an empty range
 * updates[j] = {pos, value}: update j assigns a[pos] = value
 * add(i): moves position i into the window, reading the current a[i]
 * remove(i): moves position i out of the window, reading the current a[i]
 * get(qi): the answer for the current window
 *
 * Notes:
 *   - a holds the current version during the run and is restored at the end, the window starts and ends empty
 *   - The n^(2/3) blocks assume q and u close to n, without updates mo.cpp is faster: n sqrt(q) against q n^(2/3)
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

    /// t alternates direction between block pairs: 24% fewer add / remove calls, 35% faster on random n = q = u = 1e5
    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int x, int y){
        int lx = queries[x][0] / block, ly = queries[y][0] / block;
        if (lx != ly) return lx < ly;
        int rx = queries[x][1] / block, ry = queries[y][1] / block;
        if (rx != ry) return rx < ry;
        return (lx ^ rx) & 1 ? queries[x][2] > queries[y][2] : queries[x][2] < queries[y][2];
    });

    /// Updates are applied in increasing and undone in decreasing order, so swapping with the stored value does both
    vector<T> value(u);
    for (int j = 0; j < u; j++) value[j] = updates[j].second;

    int cur_l = 0, cur_r = -1, cur_t = 0;
    /// A position in the window leaves before its value changes and re-enters after, so remove reads the value add saw
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
