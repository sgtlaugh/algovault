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
 * hld.subtree(v): the single range covering the subtree of v
 * hld.lca(u, v)
 *
 * The ranges come back in no particular order, which is fine for commutative operations (sum, min, max, xor)
 *
 * Example with a segment tree over positions:
 *     for (auto [l, r] : hld.path(u, v)) res += seg.query(l, r);
 *     for (auto [l, r] : hld.path(u, v)) seg.update(l, r, delta);
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
        vector<pair<int, int>> ranges;
        while (head[u] != head[v]){
            if (depth[head[u]] < depth[head[v]]) swap(u, v);
            ranges.push_back({pos[head[u]], pos[u]});
            u = parent[head[u]];
        }

        if (depth[u] > depth[v]) swap(u, v);
        if (pos[u] + edges <= pos[v]) ranges.push_back({pos[u] + edges, pos[v]});
        return ranges;
    }

    pair<int, int> subtree(int v) const{
        return {pos[v], pos[v] + size[v] - 1};
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

    HLD single(1);
    single.build(0);
    assert(single.path(0, 0).size() == 1 && single.path(0, 0, true).empty());

    return 0;
}
