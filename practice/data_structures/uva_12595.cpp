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

    void update_tree(vector<T>& t, int p, T v){
        for (p++; p <= n; p += p & -p) t[p] += v;
    }

    T query_tree(vector<T>& t, int p){
        T res = 0;
        for (p++; p > 0; p -= p & -p) res += t[p];
        return res;
    }
};

int main(){
    int n;

    scanf("%d", &n);
    long long c[2][3];
    for (auto& p : c) scanf("%lld %lld %lld", &p[0], &p[1], &p[2]);

    n -= 2;
    vector<long long> d1(n), d2(n);
    for (int i = 0; i < n; i++){
        long long p[3];
        scanf("%lld %lld %lld", &p[0], &p[1], &p[2]);
        d1[i] = d2[i] = 0;
        for (int k = 0; k < 3; k++){
            d1[i] += (p[k] - c[0][k]) * (p[k] - c[0][k]);
            d2[i] += (p[k] - c[1][k]) * (p[k] - c[1][k]);
        }
    }

    int q;
    scanf("%d", &q);
    vector<long long> r1(q), r2(q);
    for (int i = 0; i < q; i++) scanf("%lld %lld", &r1[i], &r2[i]);

    vector<int> by_d1(n), by_d2(n), rank(n), by_r1(q);
    iota(by_d1.begin(), by_d1.end(), 0);
    iota(by_d2.begin(), by_d2.end(), 0);
    iota(by_r1.begin(), by_r1.end(), 0);
    sort(by_d1.begin(), by_d1.end(), [&](int a, int b){ return d1[a] < d1[b]; });
    sort(by_d2.begin(), by_d2.end(), [&](int a, int b){ return d2[a] < d2[b]; });
    sort(by_r1.begin(), by_r1.end(), [&](int a, int b){ return r1[a] < r1[b]; });
    vector<long long> sorted_d2(n);
    for (int i = 0; i < n; i++) rank[by_d2[i]] = i, sorted_d2[i] = d2[by_d2[i]];

    FenwickPointUpdate<int> outside(vector<int>(n, 1));
    vector<int> res(q);
    int gone = 0;
    for (int i : by_r1){
        for (; gone < n && d1[by_d1[gone]] <= r1[i] * r1[i]; gone++) outside.update(rank[by_d1[gone]], -1);
        int k = upper_bound(sorted_d2.begin(), sorted_d2.end(), r2[i] * r2[i]) - sorted_d2.begin();
        res[i] = gone + outside.prefix(k - 1);
    }

    for (int x : res) printf("%d\n", x);
    return 0;
}
