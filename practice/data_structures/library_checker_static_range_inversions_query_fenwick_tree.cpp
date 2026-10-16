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
    vector<array<T, 2>> tree;  /// Sums of v and v*a per node; one array per node, 1.3x faster than two trees at n = 1e5 and 1.5x at n = 1e6

    FenwickFull(int n = 0) : n(n), tree(n + 1, array<T, 2>{}) {}

    void update(int l, int r, T v){
        if (l > r) return;
        assert(0 <= l && r < n);
        add_suffix(l, v);
        add_suffix(r + 1, -v);
    }

    /// A difference v at a adds v * (p - a) = v*p - v*a to the sum of the first p cells
    T prefix(int p){
        assert(-1 <= p && p < n);
        p++;
        T v = 0, va = 0;
        for (int i = p; i > 0; i -= i & -i){
            v += tree[i][0];
            va += tree[i][1];
        }
        return v * p - va;
    }

    T query(int l, int r){
        if (l > r) return 0;
        assert(0 <= l && r < n);
        return prefix(r) - prefix(l - 1);
    }

    void add_suffix(int a, T v){
        for (int i = a + 1; i <= n; i += i & -i){
            tree[i][0] += v;
            tree[i][1] += v * a;
        }
    }
};

int main(){
    int n, q;
    scanf("%d %d", &n, &q);

    vector<int> a(n);
    for (auto& x: a) scanf("%d", &x);
    vector<int> vals = a;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    for (auto& x: a) x = lower_bound(vals.begin(), vals.end(), x) - vals.begin();
    int m = vals.size();

    vector<int> l(q), r(q), order(q);
    for (int i = 0; i < q; i++) scanf("%d %d", &l[i], &r[i]);

    /// Mo's order with alternating right-end direction per block, so r sweeps back instead of resetting
    int block = max(1, (int)(n / sqrt((double)q)));
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int i, int j){
        int bi = l[i] / block, bj = l[j] / block;
        if (bi != bj) return bi < bj;
        return (bi & 1) ? r[i] > r[j] : r[i] < r[j];
    });

    /// Window is a[cl..cr - 1]; the fenwick counts its values
    FenwickPointUpdate<int> fen(m);
    vector<long long> res(q);
    long long inv = 0;
    int cl = 0, cr = 0;
    for (int i: order){
        while (cr < r[i]){
            inv += fen.query(a[cr] + 1, m - 1);
            fen.update(a[cr++], 1);
        }
        while (cl > l[i]){
            inv += fen.prefix(a[--cl] - 1);
            fen.update(a[cl], 1);
        }
        while (cr > r[i]){
            fen.update(a[--cr], -1);
            inv -= fen.query(a[cr] + 1, m - 1);
        }
        while (cl < l[i]){
            fen.update(a[cl], -1);
            inv -= fen.prefix(a[cl++] - 1);
        }
        res[i] = inv;
    }

    for (auto x: res) printf("%lld\n", x);
    return 0;
}
