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

const int MOD = 1000000007;

int main(){
    const int MAXN = 100001;

    vector<long long> pw(MAXN + 1, 1), inv(MAXN + 1, 1);
    for (int i = 1; i <= MAXN; i++){
        pw[i] = pw[i - 1] * 2 % MOD;
        inv[i] = inv[i - 1] * ((MOD + 1) / 2) % MOD;
    }

    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);
        vector<vector<int>> adj(n + 1);
        for (int i = 1; i < n; i++){
            int u, v;
            scanf("%d %d", &u, &v);
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        vector<int> parent(n + 1, 0), depth(n + 1, 0), tin(n + 1), tout(n + 1), order;
        vector<int> stack = {1};
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            tin[u] = order.size();
            order.push_back(u);
            for (int v : adj[u]){
                if (v != parent[u]) parent[v] = u, depth[v] = depth[u] + 1, stack.push_back(v);
            }
        }
        for (int i = n - 1; i >= 0; i--){
            int u = order[i];
            tout[u] = max(tout[u], tin[u]);
            if (u != 1) tout[parent[u]] = max(tout[parent[u]], tout[u]);
        }

        vector<long long> own(n + 1, 0);
        FenwickPointUpdate<long long> fen(n);
        int q;
        scanf("%d", &q);
        printf("Case %d:\n", tc);
        while (q--){
            int type, u;
            scanf("%d %d", &type, &u);
            if (type == 2){
                long long x;
                scanf("%lld", &x);
                x = (x % MOD + MOD) % MOD;
                own[u] = (own[u] + x) % MOD;
                fen.update(tin[u], x * pw[depth[u]] % MOD);
            }
            else if (tout[u] == tin[u]) printf("%lld\n", own[u]);
            else printf("%lld\n", fen.query(tin[u], tout[u]) % MOD * inv[depth[u] + 1] % MOD);
        }
    }
    return 0;
}
