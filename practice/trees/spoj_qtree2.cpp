/***
 *
 * Lowest Common Ancestor
 * LCA, distance, jumps and path folds on a weighted tree, 0-based nodes
 *
 * Complexity:
 *   - LCA: O(n log n) to build, O(log n) per lca, dist, kth_ancestor and jump, O(n log n) memory
 *   - FastLCA: O(n log n) to build, O(1) per lca and dist, O(n log n) memory
 *   - PathFold: O(n log n) to build, O(log n) per query, n log n long longs of memory
 *   - path_intersection: 6 lca calls, O(1) with FastLCA, O(log n) with LCA
 *
 * LCA tree(n) or FastLCA tree(n): an empty tree on nodes 0..n-1, both share the calls below
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * tree.build(root): roots the tree at root, 0 by default, asserts n - 1 edges and every node reached exactly once
 * tree.depth(v): the number of edges from the root to v
 * tree.dist(u, v): the sum of the edge weights on the path u - v
 * tree.lca(u, v): the lowest common ancestor of u and v
 * tree.jump(u, v, k), LCA only: the k-th vertex on the path u -> v, u for k = 0, v for the last, -1 beyond v
 * tree.kth_ancestor(v, k), LCA only: the k-th ancestor of v, v itself for k = 0, -1 if k > depth(v)
 *
 * PathFold fold(tree, op): on a built LCA, op associative and commutative on long long, like max, min or gcd
 * fold.query(u, v): op over the edge weights on the path u - v, u != v
 *
 * path_intersection(tree, a, b, c, d): the common part of paths a - b and c - d as an unordered pair of endpoints
 *
 * Notes:
 *   - path_intersection takes either structure, returns {-1, -1} for no shared vertex, {x, x} for exactly x
 *   - Use FastLCA by default, its lca is 6x faster on random trees, 20x on paths; LCA for kth_ancestor, jump, PathFold
 *   - PathFold keeps a reference to tree, so build a new fold after tree.build with another root
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LCA{
    int n, lg, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<vector<int>> up;
    vector<int> level;
    vector<long long> weight_sum;

    LCA(int n) : n(n), adj(n), level(n), weight_sum(n){
        lg = 1;
        while ((1 << lg) < n) lg++;
    }

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    /// depth < 2^lg, so levels 0..lg-1 cover every jump
    void build(int root = 0){
        assert(added_edges == n - 1);
        up.assign(lg, vector<int>(n, -1));
        vector<int> order = {root};
        up[0][root] = root, level[root] = 0, weight_sum[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (auto [v, w] : adj[u]){
                if (v == up[0][u]) continue;
                assert(up[0][v] == -1);
                up[0][v] = u, level[v] = level[u] + 1, weight_sum[v] = weight_sum[u] + w;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int k = 1; k < lg; k++){
            for (int v = 0; v < n; v++) up[k][v] = up[k - 1][up[k - 1][v]];
        }
    }

    int depth(int v) const{
        return level[v];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }

    int lca(int u, int v) const{
        if (level[u] < level[v]) swap(u, v);
        u = kth_ancestor(u, level[u] - level[v]);
        if (u == v) return u;

        /// u != v at equal depth, so level[u] >= 1, and jumps longer than level[u] reach the root from both
        for (int k = __lg(level[u]); k >= 0; k--){
            if (up[k][u] != up[k][v]) u = up[k][u], v = up[k][v];
        }

        return up[0][u];
    }

    int jump(int u, int v, int k) const{
        int l = lca(u, v), up_len = level[u] - level[l], len = up_len + level[v] - level[l];
        if (k > len) return -1;
        return k <= up_len ? kth_ancestor(u, k) : kth_ancestor(v, len - k);
    }

    int kth_ancestor(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        for (int i = 0; k; i++, k >>= 1){
            if (k & 1) v = up[i][v];
        }

        return v;
    }
};

template <typename Op>
struct PathFold{
    const LCA& tree;
    Op op;
    vector<vector<long long>> table;

    /// table[k][v] folds the 2^k edges above v; jumps past the root fold in the root's 0 and are never queried
    PathFold(const LCA& tree, Op op) : tree(tree), op(op), table(tree.lg, vector<long long>(tree.n)){
        for (int v = 0; v < tree.n; v++) table[0][v] = tree.weight_sum[v] - tree.weight_sum[tree.up[0][v]];
        for (int k = 1; k < tree.lg; k++){
            for (int v = 0; v < tree.n; v++) table[k][v] = op(table[k - 1][v], table[k - 1][tree.up[k - 1][v]]);
        }
    }

    long long query(int u, int v) const{
        assert(u != v);
        int l = tree.lca(u, v), up_u = tree.level[u] - tree.level[l], up_v = tree.level[v] - tree.level[l];
        if (up_u == 0) return climb(v, up_v);
        if (up_v == 0) return climb(u, up_u);
        return op(climb(u, up_u), climb(v, up_v));
    }

    /// op over the k >= 1 edges above v
    long long climb(int v, int k) const{
        int i = __builtin_ctz(k);
        long long res = table[i][v];

        v = tree.up[i][v];
        for (k &= k - 1; k; k &= k - 1){
            i = __builtin_ctz(k);
            res = op(res, table[i][v]);
            v = tree.up[i][v];
        }

        return res;
    }
};

/// BEGIN SHARED fast_lca
struct FastLCA{
    int n, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<int> tin, order, level;
    vector<long long> weight_sum;
    vector<vector<int>> table;

    FastLCA(int n) : n(n), adj(n), tin(n), order(n), level(n), weight_sum(n) {}

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    /// table[k][i] is the min tin over the parents of the vertices at preorder positions i..i+2^k-1
    void build(int root = 0){
        assert(added_edges == n - 1);
        vector<int> parent(n, -1);
        vector<pair<int, int>> stack = {{root, 0}};
        parent[root] = root, level[root] = 0, weight_sum[root] = 0;
        table.assign(1, vector<int>(n));
        int visited = 0;

        while (!stack.empty()){
            auto [u, parent_tin] = stack.back();
            stack.pop_back();
            tin[u] = visited, order[visited] = u, table[0][visited++] = parent_tin;
            for (auto [v, w] : adj[u]){
                if (v == parent[u]) continue;
                assert(parent[v] == -1);
                parent[v] = u, level[v] = level[u] + 1, weight_sum[v] = weight_sum[u] + w;
                stack.push_back({v, tin[u]});
            }
        }
        assert(visited == n);

        for (int k = 1; (1 << k) <= n; k++){
            table.push_back(vector<int>(n - (1 << k) + 1));
            for (int i = 0; i + (1 << k) <= n; i++){
                table[k][i] = min(table[k - 1][i], table[k - 1][i + (1 << (k - 1))]);
            }
        }
    }

    int depth(int v) const{
        return level[v];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }

    /// For tin[u] < tin[v], the shallowest vertex at preorder positions tin[u] + 1..tin[v] is a child of the lca
    int lca(int u, int v) const{
        if (u == v) return u;
        int l = tin[u], r = tin[v];
        if (l > r) swap(l, r);
        int k = __lg(r - l);
        return order[min(table[k][l + 1], table[k][r - (1 << k) + 1])];
    }
};
/// END SHARED fast_lca

/// The two deepest cross lcas are the endpoints, the paths meet iff the deepest is not above the deeper path top
template <typename Tree>
pair<int, int> path_intersection(const Tree& tree, int a, int b, int c, int d){
    int top = tree.lca(a, b), other_top = tree.lca(c, d);
    if (tree.depth(other_top) > tree.depth(top)) top = other_top;

    array<int, 4> cand = {tree.lca(a, c), tree.lca(a, d), tree.lca(b, c), tree.lca(b, d)};
    sort(cand.begin(), cand.end(), [&](int x, int y){ return tree.depth(x) > tree.depth(y); });
    if (tree.depth(cand[0]) < tree.depth(top)) return {-1, -1};
    return {cand[0], cand[1]};
}

char str[16];

int main(){
    int t, n, i, a, b, c, k;

    scanf("%d", &t);
    while (t--){
        scanf("%d", &n);
        LCA tree(n);
        for (i = 1; i < n; i++){
            scanf("%d %d %d", &a, &b, &c);
            tree.add_edge(a - 1, b - 1, c);
        }
        tree.build(0);

        while (scanf("%s", str) == 1 && str[1] != 'O'){
            scanf("%d %d", &a, &b);
            if (str[0] == 'D') printf("%lld\n", tree.dist(a - 1, b - 1));
            else{
                scanf("%d", &k);
                printf("%d\n", tree.jump(a - 1, b - 1, k - 1) + 1);
            }
        }
        puts("");
    }

    return 0;
}
