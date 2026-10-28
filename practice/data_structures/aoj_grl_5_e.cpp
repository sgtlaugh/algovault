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

int n;
vector<int> parent, heavy, head, pos;
FenwickFull<long long> fen;

void update_path(int v, long long w){
    for (; v != -1; v = parent[head[v]]) fen.update(max(pos[head[v]], 1), pos[v], w);
}

long long query_path(int v){
    long long res = 0;
    for (; v != -1; v = parent[head[v]]) res += fen.query(max(pos[head[v]], 1), pos[v]);
    return res;
}

int main(){
    scanf("%d", &n);

    vector<vector<int>> children(n);
    parent.assign(n, -1);
    for (int i = 0; i < n; i++){
        int k;
        scanf("%d", &k);
        children[i].resize(k);
        for (int& c: children[i]){
            scanf("%d", &c);
            parent[c] = i;
        }
    }

    vector<int> order = {0}, sz(n, 1);
    for (int i = 0; i < (int)order.size(); i++){
        for (int c: children[order[i]]) order.push_back(c);
    }
    heavy.assign(n, -1);
    for (int i = n - 1; i > 0; i--){
        int v = order[i], p = parent[v];
        sz[p] += sz[v];
        if (heavy[p] == -1 || sz[v] > sz[heavy[p]]) heavy[p] = v;
    }

    head.assign(n, 0);
    pos.assign(n, 0);
    int timer = 0;
    for (int v: order){
        if (v != 0 && heavy[parent[v]] == v) continue;
        for (int u = v; u != -1; u = heavy[u]){
            head[u] = v;
            pos[u] = timer++;
        }
    }

    fen = FenwickFull<long long>(n);
    int q;
    scanf("%d", &q);
    while (q--){
        int t, u;
        scanf("%d %d", &t, &u);
        if (t == 0){
            int w;
            scanf("%d", &w);
            update_path(u, w);
        }
        else printf("%lld\n", query_path(u));
    }

    return 0;
}
