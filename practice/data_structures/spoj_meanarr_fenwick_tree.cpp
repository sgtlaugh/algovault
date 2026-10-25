/***
 *
 * Fenwick Tree (Binary Indexed Tree)
 * Sums over an array under updates: point add and range sum, range add and point value, or range add and range sum
 *
 * Complexity:
 *   - build: O(n)
 *   - update, prefix, query, lower_bound: O(log n)
 *
 * FenwickPointUpdate<T> fen(n): point add, range sum, over n zeros
 * FenwickPointUpdate<T> fen(a): the same over a copy of a
 * fen.update(p, v): add v to a[p]
 * fen.prefix(p): sum of a[0..p], 0 for p = -1
 * fen.query(l, r): sum of a[l..r], 0 when l > r
 * fen.lower_bound(s): smallest p with prefix(p) >= s, n if none, every a[i] must be non-negative
 *
 * FenwickRangeUpdate<T> fen(n): range add, point value, over n zeros
 * fen.update(l, r, v): add v to a[l..r], nothing when l > r
 * fen.query(p): a[p]
 *
 * FenwickFull<T> fen(n): range add, range sum, over n zeros
 * fen.update(l, r, v): add v to a[l..r], nothing when l > r
 * fen.prefix(p): sum of a[0..p], 0 for p = -1
 * fen.query(l, r): sum of a[l..r], 0 when l > r
 *
 * Positions are 0-based, a position outside [0, n - 1] asserts
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct FenwickPointUpdate{
    int n;
    vector<T> tree;

    FenwickPointUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    /// O(n): each node pushes its finished sum into its parent
    FenwickPointUpdate(const vector<T>& a) : n(a.size()), tree(n + 1, 0){
        for (int i = 1; i <= n; i++){
            tree[i] += a[i - 1];
            if (i + (i & -i) <= n) tree[i + (i & -i)] += tree[i];
        }
    }

    void update(int p, T v){
        assert(0 <= p && p < n);
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }

    T prefix(int p){
        assert(-1 <= p && p < n);
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

    /// Binary lifting: p stays a multiple of 2 * step, so tree[p + step] is the sum of a[p .. p + step - 1]
    int lower_bound(T s){
        int p = 0, step = 1;
        while (step <= n / 2) step *= 2;

        for (; step > 0; step >>= 1){
            if (p + step <= n && tree[p + step] < s){
                p += step;
                s -= tree[p];
            }
        }

        return p;
    }
};

template <typename T>
struct FenwickRangeUpdate{
    int n;
    vector<T> tree;

    FenwickRangeUpdate(int n = 0) : n(n), tree(n + 1, 0) {}

    void update(int l, int r, T v){
        if (l > r) return;
        assert(0 <= l && r < n);
        add_suffix(l, v);
        add_suffix(r + 1, -v);
    }

    T query(int p){
        assert(0 <= p && p < n);
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += tree[p];
        return res;
    }

private:
    void add_suffix(int p, T v){
        for (p++; p <= n; p += p & -p) tree[p] += v;
    }
};

template <typename T>
struct FenwickFull{
    int n;
    vector<T> tree, aux;

    FenwickFull(int n = 0) : n(n), tree(n + 1, 0), aux(n + 1, 0) {}

    void update(int l, int r, T v){
        if (l > r) return;
        assert(0 <= l && r < n);
        update_tree(tree, l, v);
        update_tree(tree, r + 1, -v);
        update_tree(aux, l, v * l);
        update_tree(aux, r + 1, -v * (r + 1));
    }

    /// With d the difference array, sum of a[0..p] = (p + 1) * sum d[i] - sum d[i] * i; tree holds d[i], aux holds d[i] * i
    T prefix(int p){
        assert(-1 <= p && p < n);
        return query_tree(tree, p) * (p + 1) - query_tree(aux, p);
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

private:
    void update_tree(vector<T>& t, int p, T v){
        for (p++; p <= n; p += p & -p) t[p] += v;
    }

    T query_tree(vector<T>& t, int p){
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += t[p];
        return res;
    }
};

/***
 *
 * Coordinate Compression
 * Replaces each value of a vector in place with a small id, and returns the distinct values
 *
 * Complexity: O(n log n), no hashing
 *
 * compress(v): ids follow value order, returns the distinct values sorted, so vals[v[i]] is the original v[i]
 * compress(v, false): ids follow first appearance, returns the distinct values in that order
 *
 * Ids run from 0 to the number of distinct values minus 1
 *
***/



template <class T>
vector<T> compress(vector<T>& v, bool make_sorted=true){
    int n = v.size();
    vector<pair<T, int>> order(n);
    for (int i = 0; i < n; i++) order[i] = {v[i], i};
    sort(order.begin(), order.end());

    vector<T> vals;
    vector<int> ids(n);
    for (int i = 0; i < n; i++){
        if (!i || order[i - 1].first < order[i].first) vals.push_back(order[i].first);
        ids[order[i].second] = vals.size() - 1;
    }

    /// Relabel by scanning in input order instead of hashing, std::hash is the identity on ints and hackable
    if (!make_sorted){
        vector<int> relabel(vals.size(), -1);
        vector<T> first_seen;
        for (int& id: ids){
            if (relabel[id] == -1){
                relabel[id] = first_seen.size();
                first_seen.push_back(vals[id]);
            }
            id = relabel[id];
        }
        vals.swap(first_seen);
    }

    for (int i = 0; i < n; i++) v[i] = ids[i];
    return vals;
}

int main(){
    int n;
    long long k;

    scanf("%d %lld", &n, &k);
    vector<long long> p(n + 1, 0);
    for (int i = 1; i <= n; i++){
        long long a;
        scanf("%lld", &a);
        p[i] = p[i - 1] + a - k;
    }

    compress(p);
    FenwickPointUpdate<int> seen(n + 1);
    long long res = 0;
    for (int j = 0; j <= n; j++){
        res += seen.prefix(p[j]);
        seen.update(p[j], 1);
    }

    printf("%lld\n", res);
    return 0;
}
