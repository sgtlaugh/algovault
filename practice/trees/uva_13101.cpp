/***
 *
 * Heavy Light Decomposition
 * Maps every tree path to O(log n) contiguous ranges of positions, so any range structure answers path queries
 *
 * Complexity: O(n) to build, O(log n) ranges per path
 *
 * HLD hld(n); hld.add_edge(u, v); hld.build(root);
 * The graph must be a tree, build asserts that every node is reached exactly once
 * hld.pos[v]: position of v in [0, n - 1], 0-indexed like segment_tree.cpp and fenwick_tree.cpp, node values go there
 * hld.path(u, v): ranges [l, r] covering the path u..v, inclusive
 * hld.path(u, v, true): same but excludes the lca, for values on edges
 *     store the weight of edge (parent[v], v) at pos[v], the root's position stays unused
 * hld.edge_pos(u, v): position holding the value of tree edge (u, v), in either order
 * hld.ordered_path(u, v): the same ranges as {from, to} pairs listed in path order from u to v,
 *     each walked from position from to position to: from >= to on the u side (positions decrease while climbing to the lca),
 *     from <= to on the v side, a single-node range has from == to
 *     needed for non-commutative folds such as composing affine functions along the path
 * hld.ordered_path(u, v, true): same but excludes the lca
 * hld.subtree(v): the single range covering the subtree of v
 * hld.query_path(seg, u, v, edges = false): fold of the path through seg.combine, commutative combines only, T() for an empty path
 * hld.update_path(seg, u, v, x, edges = false): seg.update(l, r, x) on every range of the path
 *     seg is any size n tree with 0-indexed query(l, r), update(l, r, x) and combine(a, b), like segment_tree.cpp
 * hld.lca(u, v)
 *
 * path() returns ranges in no particular order, which is fine for commutative operations (sum, min, max, xor)
 *
 * Example with segment_tree.cpp, for another fold (max, say) pass an Ops struct:
 *     SegmentTree<long long> seg(n);
 *     for (int v = 0; v < n; v++) seg.update(hld.pos[v], hld.pos[v], value[v]);
 *     hld.update_path(seg, u, v, delta), res = hld.query_path(seg, u, v);
 *     int p = hld.edge_pos(u, v); seg.update(p, p, w), res = hld.query_path(seg, a, b, true);
 *
 * Example of a path composite, the tree keeps both the left-to-right and right-to-left fold of every node:
 *     for (auto [a, b] : hld.ordered_path(u, v)){
 *         res = a <= b ? compose(res, seg.query(a, b)) : compose(res, seg.query_reversed(b, a));
 *     }
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct HLD{
    int n;
    vector<vector<int>> adj;
    vector<int> parent, depth, heavy, head, pos, size;

    HLD(int n) : n(n), adj(n), parent(n, -1), depth(n), heavy(n, -1), head(n), pos(n), size(n, 1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void build(int root = 0){
        size.assign(n, 1), heavy.assign(n, -1), parent.assign(n, -1);  /// build may run again, e.g. with a new root
        vector<int> order = {root};
        parent[root] = -1, depth[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int v : adj[u]){
                if (v == parent[u]) continue;
                assert(v != root && parent[v] == -1);  /// a revisit means a cycle, which would grow order forever
                parent[v] = u, depth[v] = depth[u] + 1;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i > 0; i--){
            int v = order[i], p = parent[v];
            size[p] += size[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v;
        }

        /// Heavy child is pushed last so it is popped right after its parent, keeping chains and subtrees contiguous
        vector<int> stack = {root};
        head[root] = root;
        for (int cur = 0; !stack.empty(); cur++){
            int u = stack.back();
            stack.pop_back();
            pos[u] = cur;
            for (int v : adj[u]){
                if (v == parent[u] || v == heavy[u]) continue;
                head[v] = v;
                stack.push_back(v);
            }
            if (heavy[u] != -1){
                head[heavy[u]] = head[u];
                stack.push_back(heavy[u]);
            }
        }
    }

    vector<pair<int, int>> path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> ranges = ordered_path(u, v, edges);
        for (auto& e : ranges){
            if (e.first > e.second) swap(e.first, e.second);
        }

        return ranges;
    }

    vector<pair<int, int>> ordered_path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> up, down;
        while (head[u] != head[v]){
            if (depth[head[u]] >= depth[head[v]]){
                up.push_back({pos[u], pos[head[u]]});
                u = parent[head[u]];
            }
            else{
                down.push_back({pos[head[v]], pos[v]});
                v = parent[head[v]];
            }
        }

        /// u and v now share a chain, the shallower one (smaller position) is the lca
        if (pos[u] >= pos[v]){
            if (pos[u] >= pos[v] + edges) up.push_back({pos[u], pos[v] + edges});
        }
        else down.push_back({pos[u] + edges, pos[v]});

        up.insert(up.end(), down.rbegin(), down.rend());
        return up;
    }

    pair<int, int> subtree(int v) const{
        return {pos[v], pos[v] + size[v] - 1};
    }

    int edge_pos(int u, int v) const{
        return pos[depth[u] > depth[v] ? u : v];  /// edge values live on the child
    }

    /// Folds from the first range, so seg needs no identity; only u == v with edges gives an empty path
    template <typename Tree>
    auto query_path(Tree& seg, int u, int v, bool edges = false) const -> decltype(seg.query(0, 0)){
        using T = decltype(seg.query(0, 0));
        auto ranges = path(u, v, edges);
        if (ranges.empty()) return T();

        T res = seg.query(ranges[0].first, ranges[0].second);
        for (int i = 1; i < (int)ranges.size(); i++) res = seg.combine(res, seg.query(ranges[i].first, ranges[i].second));
        return res;
    }

    template <typename Tree, typename Value>
    void update_path(Tree& seg, int u, int v, Value x, bool edges = false) const{
        for (auto e : path(u, v, edges)) seg.update(e.first, e.second, x);
    }

    int lca(int u, int v) const{
        while (head[u] != head[v]){
            if (depth[head[u]] < depth[head[v]]) swap(u, v);
            u = parent[head[u]];
        }

        return depth[u] < depth[v] ? u : v;
    }
};

struct GcdTree{
    int n;
    vector<int> t;

    GcdTree(int n) : n(n), t(2 * n, 0) {}

    int combine(int a, int b) const{
        return __gcd(a, b);
    }

    void set(int p, int x){
        for (t[p += n] = x; p > 1; p >>= 1) t[p >> 1] = __gcd(t[p], t[p ^ 1]);
    }

    int query(int l, int r) const{
        int res = 0;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1){
            if (l & 1) res = __gcd(res, t[l++]);
            if (r & 1) res = __gcd(res, t[--r]);
        }

        return res;
    }
};

int main(){
    int n;

    while (scanf("%d", &n) == 1){
        vector<int> value(n);
        for (auto& x : value) scanf("%d", &x);

        HLD hld(n);
        for (int i = 1; i < n; i++){
            int u, v;
            scanf("%d %d", &u, &v);
            hld.add_edge(u, v);
        }
        hld.build(0);

        GcdTree seg(n);
        for (int v = 0; v < n; v++) seg.set(hld.pos[v], value[v]);

        int q;
        scanf("%d", &q);
        while (q--){
            int type, u, v;
            scanf("%d %d %d", &type, &u, &v);
            if (type == 2) seg.set(hld.pos[u], v);
            else printf("%d\n", hld.query_path(seg, u, v));
        }
    }

    return 0;
}
