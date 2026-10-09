/***
 *
 * Heavy Light Decomposition
 * Maps every tree path to O(log n) contiguous ranges of positions, so any range structure answers path queries
 *
 * Complexity: O(n) to build, O(log n) ranges per path
 *
 * HLD hld(n); hld.add_edge(u, v); hld.build(root);
 * hld.pos[v]: position of v in [1, n], 1-indexed like segment_tree.cpp and fenwick_tree.cpp, store node values there
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
 * hld.query_path(seg, u, v, edges = false): fold of the path through seg's own merge and t_id, commutative merges only
 * hld.update_path(seg, u, v, x, edges = false): seg.update(l, r, x) on every range of the path
 *     seg is any tree with 1-indexed query(l, r), update(l, r, x), merge(a, b) and t_id, such as segment_tree.cpp
 * hld.lca(u, v)
 *
 * path() returns ranges in no particular order, which is fine for commutative operations (sum, min, max, xor)
 *
 * Example with segment_tree.cpp, customizing its merge (sum to max, say) is the only change needed:
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
        fill(size.begin(), size.end(), 1), fill(heavy.begin(), heavy.end(), -1);  /// build may run again, e.g. with a new root
        vector<int> order = {root};
        parent[root] = -1, depth[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int v : adj[u]){
                if (v == parent[u]) continue;
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
        for (int cur = 1; !stack.empty(); cur++){
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
        for (auto& [a, b] : ranges){
            if (a > b) swap(a, b);
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

    template <typename Tree>
    auto query_path(Tree& seg, int u, int v, bool edges = false) const{
        auto res = seg.t_id;
        for (auto [l, r] : path(u, v, edges)) res = seg.merge(res, seg.query(l, r));
        return res;
    }

    template <typename Tree, typename Value>
    void update_path(Tree& seg, int u, int v, Value x, bool edges = false) const{
        for (auto [l, r] : path(u, v, edges)) seg.update(l, r, x);
    }

    int lca(int u, int v) const{
        while (head[u] != head[v]){
            if (depth[head[u]] < depth[head[v]]) swap(u, v);
            u = parent[head[u]];
        }
        return depth[u] < depth[v] ? u : v;
    }
};

int main(){
    /***
     *          0
     *        /   \
     *       1     2
     *      / \     \
     *     3   4     5
     *    /
     *   6
    ***/
    HLD hld(7);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {3, 6}}) hld.add_edge(u, v);
    hld.build(0);
    assert(hld.lca(6, 4) == 1);
    assert(hld.lca(6, 5) == 0);

    vector<long long> node_value = {1, 2, 4, 8, 16, 32, 64}, at(8);
    for (int v = 0; v < 7; v++) at[hld.pos[v]] = node_value[v];

    /// Stands in for segment_tree.cpp, with max as the merge
    struct MaxAddArray{
        vector<long long> a;
        long long t_id = LLONG_MIN;
        long long merge(long long x, long long y){ return max(x, y); }
        long long query(int l, int r){ return *max_element(a.begin() + l, a.begin() + r + 1); }
        void update(int l, int r, long long x){ for (int i = l; i <= r; i++) a[i] += x; }
    } seg{at};

    assert(hld.query_path(seg, 1, 2) == 4);     /// path 1 0 2
    hld.update_path(seg, 6, 4, 100, true);      /// edge mode skips the lca 1, adding to nodes 6 3 4
    assert(hld.query_path(seg, 6, 4) == 164);
    assert(hld.query_path(seg, 1, 1) == 2);

    auto [l, r] = hld.subtree(1);
    assert(r - l + 1 == 4);                     /// nodes 1 3 4 6
    assert(hld.edge_pos(3, 1) == hld.pos[3]);   /// an edge value lives on its child

    /// ordered_path walks positions in path order, needed for non-commutative folds
    vector<int> node_at(8), walk;
    for (int v = 0; v < 7; v++) node_at[hld.pos[v]] = v;
    for (auto [a, b] : hld.ordered_path(5, 6)){
        int step = a <= b ? 1 : -1;
        for (int i = a; i != b + step; i += step) walk.push_back(node_at[i]);
    }
    assert((walk == vector<int>{5, 2, 0, 1, 3, 6}));
    return 0;
}
