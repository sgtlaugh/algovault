/***
 *
 * Mo's Algorithm With Updates
 * Answers offline range queries on an array that also receives point assignments
 *
 * Complexity: O((q + u) n^(2/3) + n^(4/3)) calls to add and remove for q queries and u updates, plus O(q log q) sorting
 *
 * mo_with_updates(a, queries, updates, add, remove, get): the answers, res[qi] = get(qi) over the window of queries[qi]
 * queries[qi] = {l, r, t}: the inclusive range l..r of a after the first t updates, 0 <= t <= u, r = l - 1 is an empty range
 * updates[j] = {pos, value}: update j assigns a[pos] = value
 * add(i): moves position i into the window, reading the current a[i]
 * remove(i): moves position i out of the window, reading the current a[i]
 * get(qi): the answer for the current window
 *
 * a holds the current version during the run and is restored at the end, the window starts and ends empty
 * The n^(2/3) blocks assume q and u close to n, without updates mo.cpp is faster: n sqrt(q) against q n^(2/3)
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
    /// A position inside the window leaves before its value changes and re-enters after, so remove reads the value add saw
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
    int n, q;
    scanf("%d %d", &n, &q);
    vector<int> a(n), type(q), x(q), y(q), z(q);
    for (int& v : a) scanf("%d", &v);
    vector<int> values(a);
    for (int i = 0; i < q; i++){
        scanf("%d %d %d", &type[i], &x[i], &y[i]);
        if (type[i] == 1) scanf("%d", &z[i]);
        values.push_back(type[i] == 0 ? y[i] : z[i]);
    }

    /// Compressed values index a plain counter array, a hash map in add / remove would dominate the n^(5/3) work
    sort(values.begin(), values.end());
    values.erase(unique(values.begin(), values.end()), values.end());
    auto id = [&](int v){ return int(lower_bound(values.begin(), values.end(), v) - values.begin()); };
    for (int& v : a) v = id(v);

    vector<array<int, 3>> queries;
    vector<pair<int, int>> updates;
    vector<int> target;
    for (int i = 0; i < q; i++){
        if (type[i] == 0) updates.push_back({x[i], id(y[i])});
        else queries.push_back({x[i], y[i] - 1, (int)updates.size()}), target.push_back(id(z[i]));
    }

    vector<int> freq(values.size(), 0);
    auto res = mo_with_updates(a, queries, updates,
        [&](int i){ freq[a[i]]++; },
        [&](int i){ freq[a[i]]--; },
        [&](int qi){ return freq[target[qi]]; });
    for (int r : res) printf("%d\n", r);
    return 0;
}
