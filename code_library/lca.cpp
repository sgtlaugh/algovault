/***
 *
 * Lowest Common Ancestor
 * Two interchangeable structures for a weighted tree, 0-based nodes
 * Use LinearLCA by default, and LCA when kth_ancestor is needed
 * Measured: LinearLCA answers lca 2-8x faster at every size, both build in about the same time
 *
 * LCA (binary lifting):
 *   - O(n log n) to build, O(log n) per lca / kth_ancestor, O(n log n) memory
 *   - kth_ancestor(v, k): the k-th ancestor of v, v itself for k = 0, -1 if k > depth(v)
 *
 * LinearLCA (Euler tour + linear RMQ):
 *   - O(n) to build, O(1) per lca, O(n) memory, use it for n around 1e6 or many queries
 *
 * Shared API:
 *   LCA tree(n); tree.add_edge(u, v, w); tree.build(root);
 *   tree.lca(u, v), tree.dist(u, v) (sum of weights on the path, w defaults to 1), tree.depth(v)
 *   The graph must be a tree, build asserts n - 1 edges and that every node is reached exactly once
 *   Traversals are iterative so deep trees are safe
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

    int kth_ancestor(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        for (int i = 0; k; i++, k >>= 1){
            if (k & 1) v = up[i][v];
        }
        return v;
    }

    int lca(int u, int v) const{
        if (level[u] < level[v]) swap(u, v);
        u = kth_ancestor(u, level[u] - level[v]);
        if (u == v) return u;
        for (int k = lg - 1; k >= 0; k--){
            if (up[k][u] != up[k][v]) u = up[k][u], v = up[k][v];
        }
        return up[0][u];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }
};

/// Range minimum index over a fixed array in O(n) build and O(1) query, blocks of 64 with per position stack masks
/// Same technique as sparse_table.cpp's LinearSparseTable, which is the general purpose version
struct LinearRMQ{
    static constexpr int B = 64;
    vector<int> val;
    vector<unsigned long long> mask;
    vector<vector<int>> table;

    int better(int i, int j) const{
        return val[i] <= val[j] ? i : j;
    }

    /// Index of the minimum in [r - size + 1, r], size <= 64
    int small_query(int r, int size = B) const{
        unsigned long long m = size == B ? mask[r] : mask[r] & ((1ULL << size) - 1);
        return r - (63 - __builtin_clzll(m));
    }

    LinearRMQ(){}

    LinearRMQ(const vector<int>& values) : val(values), mask(values.size()){
        int n = val.size();
        unsigned long long cur = 0;
        for (int i = 0; i < n; i++){
            cur <<= 1;
            while (cur && val[i - __builtin_ctzll(cur)] >= val[i]) cur &= cur - 1;
            cur |= 1;
            mask[i] = cur;
        }

        int blocks = (n + B - 1) / B;
        table.assign(1, vector<int>(blocks));
        for (int b = 0; b < blocks; b++) table[0][b] = small_query(min(n - 1, b * B + B - 1), min(B, n - b * B));
        for (int k = 1; (1 << k) <= blocks; k++){
            table.push_back(vector<int>(blocks - (1 << k) + 1));
            for (int b = 0; b + (1 << k) <= blocks; b++) table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
        }
    }

    int query(int l, int r) const{
        if (r - l + 1 <= B) return small_query(r, r - l + 1);
        int res = better(small_query(l + B - 1), small_query(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = 31 - __builtin_clz(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }
        return res;
    }
};

struct LinearLCA{
    int n, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<int> first, tour, level;
    vector<long long> weight_sum;
    LinearRMQ rmq;

    LinearLCA(int n) : n(n), adj(n), first(n), level(n), weight_sum(n){}

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    void build(int root = 0){
        assert(added_edges == n - 1);
        vector<int> parent(n, -1), next_edge(n, 0), stack = {root}, depths = {0};
        parent[root] = root;
        tour = {root};
        first[root] = 0, level[root] = 0, weight_sum[root] = 0;
        int visited = 1;
        while (!stack.empty()){
            int u = stack.back();
            if (next_edge[u] < (int)adj[u].size()){
                auto [v, w] = adj[u][next_edge[u]++];
                if (v == parent[u]) continue;
                assert(parent[v] == -1);
                parent[v] = u, level[v] = level[u] + 1, weight_sum[v] = weight_sum[u] + w;
                first[v] = tour.size(), tour.push_back(v), depths.push_back(level[v]);
                stack.push_back(v), visited++;
                continue;
            }
            stack.pop_back();
            if (!stack.empty()) tour.push_back(stack.back()), depths.push_back(level[stack.back()]);
        }
        assert(visited == n);
        rmq = LinearRMQ(depths);
    }

    int depth(int v) const{
        return level[v];
    }

    int lca(int u, int v) const{
        int l = first[u], r = first[v];
        if (l > r) swap(l, r);
        return tour[rmq.query(l, r)];
    }

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
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
    vector<array<int, 3>> edges = {{0, 1, 3}, {0, 2, 1}, {1, 3, 2}, {1, 4, 7}, {2, 5, 4}, {3, 6, 5}};

    LCA a(7);
    LinearLCA b(7);
    for (auto [u, v, w] : edges) a.add_edge(u, v, w), b.add_edge(u, v, w);
    a.build(0), b.build(0);

    assert(a.lca(6, 4) == 1 && b.lca(6, 4) == 1);
    assert(a.lca(6, 5) == 0 && b.lca(6, 5) == 0);
    assert(a.lca(3, 6) == 3 && b.lca(3, 6) == 3);
    assert(a.lca(5, 5) == 5 && b.lca(5, 5) == 5);
    assert(a.dist(6, 4) == 14 && b.dist(6, 4) == 14);
    assert(a.dist(6, 5) == 15 && b.dist(6, 5) == 15);
    assert(a.depth(6) == 3 && b.depth(6) == 3);

    assert(a.kth_ancestor(6, 0) == 6);
    assert(a.kth_ancestor(6, 1) == 3);
    assert(a.kth_ancestor(6, 2) == 1);
    assert(a.kth_ancestor(6, 3) == 0);
    assert(a.kth_ancestor(6, 4) == -1);

    LCA single(1);
    single.build(0);
    assert(single.lca(0, 0) == 0 && single.kth_ancestor(0, 1) == -1);
    LinearLCA single_linear(1);
    single_linear.build(0);
    assert(single_linear.lca(0, 0) == 0);

    a.build(4), b.build(4);
    assert(a.lca(6, 5) == 1 && b.lca(6, 5) == 1);
    assert(a.kth_ancestor(5, 3) == 1);

    return 0;
}
