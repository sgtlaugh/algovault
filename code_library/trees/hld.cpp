/***
 *
 * Heavy Light Decomposition
 * Maps every tree path to O(log n) contiguous ranges of positions, so any range structure answers path queries
 *
 * Complexity: O(n) to build, O(log n) ranges per path
 *
 * HLD hld(n); hld.add_edge(u, v); hld.build(root);
 * hld.pos[v]: position of v in [0, n), store node values at these positions in your segment tree / fenwick
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
 * hld.lca(u, v)
 *
 * path() returns ranges in no particular order, which is fine for commutative operations (sum, min, max, xor)
 *
 * Example with a segment tree over positions:
 *     for (auto [l, r] : hld.path(u, v)) res += seg.query(l, r);
 *     for (auto [l, r] : hld.path(u, v)) seg.update(l, r, delta);
 *     seg.update(hld.edge_pos(u, v), w), then for (auto [l, r] : hld.path(a, b, true)) res += seg.query(l, r);
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

    vector<long long> node_value = {1, 2, 4, 8, 16, 32, 64}, at(7);
    for (int v = 0; v < 7; v++) at[hld.pos[v]] = node_value[v];
    auto path_sum = [&](int u, int v, bool edges){
        long long s = 0;
        for (auto [l, r] : hld.path(u, v, edges)){
            for (int i = l; i <= r; i++) s += at[i];
        }
        return s;
    };

    assert(path_sum(6, 4, false) == 64 + 8 + 2 + 16);
    assert(path_sum(6, 5, false) == 64 + 8 + 2 + 1 + 4 + 32);
    assert(path_sum(3, 3, false) == 8);
    assert(path_sum(6, 4, true) == 64 + 8 + 16);
    assert(path_sum(3, 3, true) == 0);
    assert(path_sum(6, 5, true) == 64 + 8 + 2 + 4 + 32);

    auto [l, r] = hld.subtree(1);
    long long s = 0;
    for (int i = l; i <= r; i++) s += at[i];
    assert(s == 2 + 8 + 16 + 64);
    assert(hld.subtree(5).first == hld.subtree(5).second);

    assert(hld.lca(6, 4) == 1);
    assert(hld.lca(6, 5) == 0);
    assert(hld.lca(3, 6) == 3);

    assert(hld.edge_pos(1, 3) == hld.pos[3] && hld.edge_pos(3, 1) == hld.pos[3] && hld.edge_pos(0, 2) == hld.pos[2]);

    vector<int> node_at(7);
    for (int v = 0; v < 7; v++) node_at[hld.pos[v]] = v;
    auto ordered_nodes = [&](int u, int v, bool edges){
        vector<int> nodes;
        for (auto [a, b] : hld.ordered_path(u, v, edges)){
            int step = a <= b ? 1 : -1;
            for (int i = a; i != b + step; i += step) nodes.push_back(node_at[i]);
        }
        return nodes;
    };

    assert((ordered_nodes(6, 4, false) == vector<int>{6, 3, 1, 4}));
    assert((ordered_nodes(4, 6, false) == vector<int>{4, 1, 3, 6}));
    assert((ordered_nodes(6, 5, false) == vector<int>{6, 3, 1, 0, 2, 5}));
    assert((ordered_nodes(5, 6, false) == vector<int>{5, 2, 0, 1, 3, 6}));
    assert((ordered_nodes(6, 1, false) == vector<int>{6, 3, 1}));
    assert((ordered_nodes(1, 6, false) == vector<int>{1, 3, 6}));
    assert((ordered_nodes(6, 4, true) == vector<int>{6, 3, 4}));
    assert((ordered_nodes(5, 6, true) == vector<int>{5, 2, 1, 3, 6}));
    assert((ordered_nodes(1, 6, true) == vector<int>{3, 6}));
    assert((ordered_nodes(3, 3, false) == vector<int>{3}));
    assert((ordered_nodes(3, 3, true) == vector<int>{}));

    vector<pair<long long, long long>> affine = {{1, -1}, {1, 5}, {2, 0}, {3, 0}, {10, 0}, {1, 3}, {2, 1}}, func_at(7);
    for (int v = 0; v < 7; v++) func_at[hld.pos[v]] = affine[v];
    auto path_apply = [&](int u, int v, bool edges, long long x){
        for (auto [a, b] : hld.ordered_path(u, v, edges)){
            int step = a <= b ? 1 : -1;
            for (int i = a; i != b + step; i += step) x = func_at[i].first * x + func_at[i].second;
        }
        return x;
    };

    assert(path_apply(6, 4, false, 1) == 140);
    assert(path_apply(4, 6, false, 1) == 91);
    assert(path_apply(6, 4, true, 1) == 90);
    assert(path_apply(6, 5, false, 1) == 29);
    assert(path_apply(5, 6, false, 1) == 73);

    HLD single(1);
    single.build(0);
    assert(single.path(0, 0).size() == 1 && single.path(0, 0, true).empty());
    assert((single.ordered_path(0, 0) == vector<pair<int, int>>{{0, 0}}) && single.ordered_path(0, 0, true).empty());

    return 0;
}
