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

    vector<long long> a(n);
    for (auto& x: a) scanf("%lld", &x);

    vector<vector<int>> adj(n);
    for (int i = 1; i < n; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    /// Centroid c keeps one fenwick over distances to c for its whole component, and one per child piece so a query
    /// from inside that piece can subtract what the lower centroids already counted
    vector<FenwickPointUpdate<long long>> fens;
    vector<vector<array<int, 3>>> entries(n);  /// Per vertex, per centroid ancestor: {centroid fenwick, distance, piece fenwick or -1}
    vector<bool> removed(n, false);
    vector<int> par(n), sz(n), dist(n), piece(n), stk = {0};

    while (!stk.empty()){
        int root = stk.back();
        stk.pop_back();

        vector<int> comp = {root};
        par[root] = -1;
        for (size_t i = 0; i < comp.size(); i++){
            int v = comp[i];
            for (int w: adj[v]){
                if (w == par[v] || removed[w]) continue;
                par[w] = v;
                comp.push_back(w);
            }
        }
        int total = comp.size(), c = root;
        for (int i = total - 1; i >= 0; i--){
            int v = comp[i], biggest = 0;
            sz[v] = 1;
            for (int w: adj[v]){
                if (w == par[v] || removed[w]) continue;
                sz[v] += sz[w];
                biggest = max(biggest, sz[w]);
            }
            if (max(biggest, total - sz[v]) <= total / 2) c = v;
        }

        vector<int> bfs = {c};
        par[c] = -1;
        dist[c] = 0;
        piece[c] = -1;
        int pieces = 0;
        for (size_t i = 0; i < bfs.size(); i++){
            int v = bfs[i];
            for (int w: adj[v]){
                if (w == par[v] || removed[w]) continue;
                par[w] = v;
                dist[w] = dist[v] + 1;
                piece[w] = v == c ? pieces++ : piece[v];
                bfs.push_back(w);
            }
        }

        vector<long long> whole(dist[bfs.back()] + 1, 0);
        vector<vector<long long>> parts(pieces);
        for (int v: bfs){
            whole[dist[v]] += a[v];
            if (piece[v] == -1) continue;
            auto& p = parts[piece[v]];
            if ((int)p.size() <= dist[v]) p.resize(dist[v] + 1, 0);
            p[dist[v]] += a[v];
        }

        int base = fens.size();
        fens.emplace_back(whole);
        for (auto& p: parts) fens.emplace_back(p);
        for (int v: bfs) entries[v].push_back({base, dist[v], piece[v] == -1 ? -1 : base + 1 + piece[v]});

        removed[c] = true;
        for (int w: adj[c]){
            if (!removed[w]) stk.push_back(w);
        }
    }

    auto range = [&](int id, int lo, int hi){
        lo = max(lo, 0);
        hi = min(hi, fens[id].n - 1);
        return fens[id].query(lo, hi);
    };

    while (q--){
        int t, p;
        scanf("%d %d", &t, &p);
        if (t == 0){
            long long x;
            scanf("%lld", &x);
            for (auto& e: entries[p]){
                fens[e[0]].update(e[1], x);
                if (e[2] != -1) fens[e[2]].update(e[1], x);
            }
        }
        else{
            int l, r;
            scanf("%d %d", &l, &r);
            long long res = 0;
            for (auto& e: entries[p]){
                res += range(e[0], l - e[1], r - 1 - e[1]);
                if (e[2] != -1) res -= range(e[2], l - e[1], r - 1 - e[1]);
            }
            printf("%lld\n", res);
        }
    }
    return 0;
}
