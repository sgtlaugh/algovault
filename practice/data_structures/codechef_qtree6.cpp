/***
 *
 * Link-Cut Tree
 * Dynamic rooted forest with link, cut, re-rooting, LCA and path aggregates
 *
 * Complexity: amortized O(log n) per operation
 *
 * LinkCutTree<T> lct(n): nodes 0..n - 1 and no edges, values start at T() and paths fold with +
 * LinkCutTree<T, Combine> lct(n, identity, combine): values start at identity and paths fold with combine
 * lct.link(u, v): adds edge u-v with u as a child of v, false and no change if u and v are already connected
 * lct.cut(u, v): removes edge u-v, false if there is no such edge
 * lct.connected(u, v): whether u and v are in the same tree
 * lct.find_root(u): root of u's tree
 * lct.parent(u): parent of u under the current root, -1 for the root
 * lct.make_root(u): re-roots u's tree at u
 * lct.lca(u, v): lowest common ancestor under the current root, -1 if u and v are not connected
 * lct.set(u, value): sets the value of u
 * lct.query(u, v): combine of the values on the path u..v, both ends included, u and v must be connected
 *
 * combine must be associative and commutative with identity as its neutral element, since re-rooting reverses paths
 * link re-roots u's tree at u and v's tree keeps its root, so link(child, parent) builds a rooted forest for lca
 * After cut the piece holding the old root keeps it, the other piece is rooted at its endpoint of the edge
 * query leaves the root where it was
 * Everything is iterative, deep paths are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Combine = plus<T>>
struct LinkCutTree{
    struct Node{
        /// par of an auxiliary splay tree's root is the path-parent pointer into the splay tree above it
        int ch[2], par;
        bool flip;
        T value, agg;
    };

    T identity;
    Combine combine;
    vector<Node> nodes;

    LinkCutTree(int n, T identity = T(), Combine combine = Combine())
        : identity(identity), combine(combine), nodes(n, Node{{-1, -1}, -1, false, identity, identity}) {}

    int parent(int u){
        access(u);
        int x = nodes[u].ch[0];
        if (x == -1) return -1;

        for (push(x); nodes[x].ch[1] != -1; push(x)) x = nodes[x].ch[1];
        splay(x);
        return x;
    }

    int find_root(int u){
        access(u);
        int x = u;
        for (push(x); nodes[x].ch[0] != -1; push(x)) x = nodes[x].ch[0];

        splay(x);
        return x;
    }

    bool connected(int u, int v){
        return find_root(u) == find_root(v);
    }

    void make_root(int u){
        access(u);
        nodes[u].flip ^= 1;
    }

    bool link(int u, int v){
        if (connected(u, v)) return false;
        make_root(u);
        nodes[u].par = v;
        return true;
    }

    bool cut(int u, int v){
        if (parent(u) != v) swap(u, v);
        if (parent(u) != v) return false;

        access(u);
        nodes[nodes[u].ch[0]].par = -1;
        nodes[u].ch[0] = -1;
        pull(u);
        return true;
    }

    int lca(int u, int v){
        if (!connected(u, v)) return -1;
        access(u);
        return access(v);
    }

    void set(int u, T value){
        splay(u);
        nodes[u].value = value;
        pull(u);
    }

    T query(int u, int v){
        assert(connected(u, v));
        int root = find_root(u);

        make_root(u);
        access(v);
        T res = nodes[v].agg;

        make_root(root);
        return res;
    }

    bool is_root(int x) const{
        int p = nodes[x].par;
        return p == -1 || (nodes[p].ch[0] != x && nodes[p].ch[1] != x);
    }

    T agg_of(int x) const{
        return x == -1 ? identity : nodes[x].agg;
    }

    void push(int x){
        if (!nodes[x].flip) return;
        swap(nodes[x].ch[0], nodes[x].ch[1]);
        for (int c : nodes[x].ch){
            if (c != -1) nodes[c].flip ^= 1;
        }
        nodes[x].flip = false;
    }

    void pull(int x){
        nodes[x].agg = combine(combine(agg_of(nodes[x].ch[0]), nodes[x].value), agg_of(nodes[x].ch[1]));
    }

    void rotate(int x){
        int p = nodes[x].par, g = nodes[p].par, d = nodes[p].ch[1] == x, child = nodes[x].ch[!d];

        if (!is_root(p)) nodes[g].ch[nodes[g].ch[1] == p] = x;
        nodes[x].par = g;
        nodes[p].ch[d] = child;
        if (child != -1) nodes[child].par = p;
        nodes[x].ch[!d] = p;
        nodes[p].par = x;

        pull(p);
        pull(x);
    }

    /// Pushing only g, p, x is enough: a flip pending higher up reverses the whole subtree and commutes with rotations below
    void splay(int x){
        push(x);
        while (!is_root(x)){
            int p = nodes[x].par, g = nodes[p].par;
            if (!is_root(p)) push(g);
            push(p);
            push(x);
            if (!is_root(p)) rotate((nodes[g].ch[1] == p) == (nodes[p].ch[1] == x) ? p : x);
            rotate(x);
        }
    }

    /// Makes root..x the preferred path, x ends as the root of its splay tree; returns the last node joined
    int access(int x){
        int last = -1;
        for (int y = x; y != -1; y = nodes[y].par){
            splay(y);
            nodes[y].ch[1] = last;
            pull(y);
            last = y;
        }

        splay(x);
        return last;
    }
};

struct Fenwick{
    vector<long long> bit;

    Fenwick(int n) : bit(n + 1, 0) {}

    void add(int i, long long x){
        for (i++; i < (int)bit.size(); i += i & -i) bit[i] += x;
    }

    long long sum(int i) const{
        long long res = 0;
        for (; i > 0; i -= i & -i) res += bit[i];
        return res;
    }
};

int main(){
    int n;

    scanf("%d", &n);
    vector<vector<int>> adj(n + 1);
    for (int i = 1; i < n; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> par(n + 1, -1), depth(n + 1, 0), order = {1}, tin(n + 1), tout(n + 1), sz(n + 1, 1);
    par[1] = 0;
    for (int i = 0; i < (int)order.size(); i++){
        int u = order[i];
        for (int v : adj[u]){
            if (v != par[u]) par[v] = u, depth[v] = depth[u] + 1, order.push_back(v);
        }
    }
    for (int i = n - 1; i > 0; i--) sz[par[order[i]]] += sz[order[i]];

    vector<int> stk = {1};
    for (int timer = 0; !stk.empty();){
        int u = stk.back();
        stk.pop_back();
        tin[u] = timer++;
        tout[u] = tin[u] + sz[u];
        for (int v : adj[u]){
            if (v != par[u]) stk.push_back(v);
        }
    }

    const long long NONE = LLONG_MAX;
    auto shallowest = [](long long a, long long b){ return min(a, b); };
    vector<LinkCutTree<long long, decltype(shallowest)>> forest(2, LinkCutTree<long long, decltype(shallowest)>(n + 1, NONE, shallowest));
    vector<Fenwick> delta(2, Fenwick(n));
    vector<int> color(n + 1, 0);
    auto key = [&](int v){ return (long long)depth[v] << 20 | v; };
    auto cnt = [&](int c, int v){ return (c == 0 ? sz[v] : 1) + delta[c].sum(tout[v]) - delta[c].sum(tin[v]); };
    auto path_add = [&](int c, int from, int top, long long x){
        delta[c].add(tin[from], x);
        if (top > 1) delta[c].add(tin[par[top]], -x);
    };

    for (int v = 1; v <= n; v++){
        forest[0].set(v, key(v));
        forest[0].link(v, par[v]);
    }

    int m;
    scanf("%d", &m);
    while (m--){
        int t, u;
        scanf("%d %d", &t, &u);
        int c = color[u];
        if (t == 0){
            int top = forest[c].query(u, forest[c].find_root(u)) & ((1 << 20) - 1);
            printf("%lld\n", cnt(c, top));
            continue;
        }

        int p = par[u];
        long long moved = cnt(c, u);
        forest[c].cut(u, p);
        forest[c].set(u, NONE);
        if (p != 0) path_add(c, p, forest[c].find_root(p), -moved);

        c ^= 1;
        color[u] = c;
        moved = cnt(c, u);
        forest[c].set(u, key(u));
        forest[c].link(u, p);
        if (p != 0) path_add(c, p, forest[c].find_root(p), moved);
    }
    return 0;
}
