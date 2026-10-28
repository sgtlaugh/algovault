/***
 *
 * Heavy Light Decomposition
 * Maps every tree path to O(log n) contiguous ranges of positions, so any range structure answers path queries
 *
 * Complexity: O(n) to build, O(log n) ranges per path
 *
 * HLD hld(n): an empty tree on nodes 0..n-1
 * hld.add_edge(u, v): adds the edge u - v
 * hld.build(root): roots the tree at root, 0 by default, asserts that every node is reached exactly once
 * hld.update_path(seg, u, v, x, edges = false): seg.update(l, r, x) on every range of the path u - v
 * hld.edge_pos(u, v): the position holding the value of tree edge (u, v), in either order, asserts u and v are adjacent
 * hld.lca(u, v): the lowest common ancestor of u and v
 * hld.ordered_path(u, v, edges = false): the ranges of the path as {from, to} pairs, in path order from u to v
 * hld.path(u, v, edges = false): ranges [l, r] covering the path u - v, in no particular order
 * hld.query_path(seg, u, v, edges = false): the fold of the path through seg.combine, T() for an empty path
 * hld.subtree(v): the single range covering the subtree of v
 * hld.pos[v]: the position of v in [0, n - 1], where the value of node v goes
 * hld.parent[v]: the parent of v, -1 for the root
 * hld.depth[v]: the number of edges from the root to v
 *
 * Notes:
 *   - seg is any size n tree with 0-indexed query(l, r), update(l, r, x) and combine(a, b), like segment_tree.cpp
 *   - edges = true excludes the lca, for values on edges: store the weight of edge (parent[v], v) at pos[v]
 *   - ordered_path walks each range from position from to to: from >= to on the u side, from <= to on the v side
 *   - path and query_path suit commutative folds only (sum, min, max, xor), ordered_path suits any, like affine maps
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct HLD{
    struct Chain{ int head_pos, head_parent; };

    int n;
    vector<pair<int, int>> edge_list;
    vector<int> parent, depth, heavy, head, pos, size;
    vector<Chain> chain;  /// v's head pos and parent in one cache line: lca 2-3x faster on random trees

    HLD(int n) : n(n), parent(n, -1), depth(n), heavy(n, -1), head(n), pos(n), size(n, 1), chain(n) {}

    void add_edge(int u, int v){
        edge_list.push_back({u, v});
    }

    void build(int root = 0){
        size.assign(n, 1), heavy.assign(n, -1), parent.assign(n, -1);  /// build may run again, e.g. with a new root

        /// One flat array instead of a vector per node: the neighbors of u are adj[start[u]..start[u + 1] - 1]
        vector<int> start(n + 1, 0), adj(2 * edge_list.size());
        for (auto [u, v] : edge_list) start[u]++, start[v]++;
        for (int i = 1; i <= n; i++) start[i] += start[i - 1];
        for (auto [u, v] : edge_list) adj[--start[u]] = v, adj[--start[v]] = u;

        vector<int> order = {root};
        parent[root] = -1, depth[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j];
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
            chain[u] = {pos[head[u]], parent[head[u]]};  /// head[u] is u or was popped before u
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j];
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

    template <typename Tree, typename Value>
    void update_path(Tree& seg, int u, int v, Value x, bool edges = false) const{
        for (auto [l, r] : path(u, v, edges)) seg.update(l, r, x);
    }

    int edge_pos(int u, int v) const{
        assert(parent[u] == v || parent[v] == u);
        return pos[depth[u] > depth[v] ? u : v];  /// edge values live on the child
    }

    /// The head with the larger position is never an ancestor of the other node, so its chain is safe to leave
    int lca(int u, int v) const{
        while (chain[u].head_pos != chain[v].head_pos){
            if (chain[u].head_pos < chain[v].head_pos) swap(u, v);
            u = chain[u].head_parent;
        }
        return pos[u] < pos[v] ? u : v;
    }

    vector<pair<int, int>> ordered_path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> up, down;

        while (chain[u].head_pos != chain[v].head_pos){
            if (chain[u].head_pos > chain[v].head_pos){
                up.push_back({pos[u], chain[u].head_pos});
                u = chain[u].head_parent;
            }
            else{
                down.push_back({chain[v].head_pos, pos[v]});
                v = chain[v].head_parent;
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

    vector<pair<int, int>> path(int u, int v, bool edges = false) const{
        vector<pair<int, int>> ranges = ordered_path(u, v, edges);
        for (auto& [a, b] : ranges){
            if (a > b) swap(a, b);
        }
        return ranges;
    }

    /// Folds from the first range, so seg needs no identity; only u == v with edges gives an empty path
    template <typename Tree>
    auto query_path(Tree& seg, int u, int v, bool edges = false) const{
        using T = decltype(seg.query(0, 0));
        auto ranges = path(u, v, edges);
        if (ranges.empty()) return T();

        T res = seg.query(ranges[0].first, ranges[0].second);
        for (int i = 1; i < (int)ranges.size(); i++){
            res = seg.combine(res, seg.query(ranges[i].first, ranges[i].second));
        }
        return res;
    }

    pair<int, int> subtree(int v) const{
        return {pos[v], pos[v] + size[v] - 1};
    }
};

#define MAX 10010

namespace seg{
    int n, tree[2 * MAX];

    void init(int m){
        n = m;
        for (int i = 0; i < 2 * n; i++) tree[i] = INT_MIN;
    }

    void update(int i, int x){
        for (tree[i += n] = x; i > 1; i >>= 1) tree[i >> 1] = max(tree[i], tree[i ^ 1]);
    }

    int query(int l, int r){
        int res = INT_MIN;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1){
            if (l & 1) res = max(res, tree[l++]);
            if (r & 1) res = max(res, tree[--r]);
        }
        return res;
    }
}

int n, eu[MAX], ev[MAX], ec[MAX];
char str[20];

int main(){
    int t, i, a, b, res;

    scanf("%d", &t);
    while (t--){
        scanf("%d", &n);
        HLD hld(n);
        for (i = 1; i < n; i++){
            scanf("%d %d %d", &eu[i], &ev[i], &ec[i]);
            hld.add_edge(--eu[i], --ev[i]);
        }
        hld.build(0);

        seg::init(n);
        for (i = 1; i < n; i++) seg::update(hld.edge_pos(eu[i], ev[i]), ec[i]);

        while (scanf("%s", str) == 1 && str[0] != 'D'){
            scanf("%d %d", &a, &b);
            if (str[0] == 'C') seg::update(hld.edge_pos(eu[a], ev[a]), b);
            else{
                res = INT_MIN;
                for (auto p : hld.path(a - 1, b - 1, true)) res = max(res, seg::query(p.first, p.second));
                printf("%d\n", a == b ? 0 : res);
            }
        }
    }
    return 0;
}
