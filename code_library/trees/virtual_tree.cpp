/***
 *
 * Virtual Tree (Auxiliary Tree)
 * Compresses a rooted tree to a vertex subset plus the pairwise LCAs needed to keep its shape, 0-based nodes
 *
 * Complexity: O(n) to build, O(k log k) per compress of k vertices
 *
 * VirtualTree tree(n): an empty tree on nodes 0..n-1
 * tree.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * tree.build(root): roots the tree at root, 0 by default, asserts n - 1 edges and every node reached exactly once
 * tree.compress(nodes): the virtual tree of nodes as a list of VirtualTree::Node {v, parent, len}
 * tree.depth(v): the number of edges from the root to v
 * tree.dist(u, v): the sum of the edge weights on the path u - v
 * tree.lca(u, v): the lowest common ancestor of u and v
 *
 * Notes:
 *   - parent is the index of v's parent in the list, -1 for the root at index 0
 *   - len is the original length of the compressed edge, dist(v, list[parent].v), 0 for the root
 *   - The list is in DFS order, so parents come before children and every subtree is a contiguous range
 *   - The list has at most 2k - 1 vertices, duplicates in nodes are ignored, an empty nodes gives an empty list
 *   - Traversals are iterative so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY linear_lca from code_library/trees/lca.cpp
/// Range minimum index in O(n) build and O(1) query, the technique of sparse_table.cpp's LinearSparseTable
struct LinearRMQ{
    static constexpr int B = 64;
    vector<int> val;
    vector<unsigned long long> mask;
    vector<vector<int>> table;

    LinearRMQ() {}

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
            for (int b = 0; b + (1 << k) <= blocks; b++){
                table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
            }
        }
    }

    int query(int l, int r) const{
        if (r - l + 1 <= B) return small_query(r, r - l + 1);
        int res = better(small_query(l + B - 1), small_query(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = __lg(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }
        return res;
    }

    /// Index of the minimum in [r - size + 1, r], size <= 64
    int small_query(int r, int size = B) const{
        unsigned long long m = size == B ? mask[r] : mask[r] & ((1ULL << size) - 1);
        /// Not __lg: GCC 14's __lg tests for 0, which costs LinearLCA.lca 4-8%
        return r - (63 - __builtin_clzll(m));
    }

    int better(int i, int j) const{
        return val[i] <= val[j] ? i : j;
    }
};

struct LinearLCA{
    int n, added_edges = 0;
    vector<vector<pair<int, long long>>> adj;
    vector<int> first, tour, level;
    vector<long long> weight_sum;
    LinearRMQ rmq;

    LinearLCA(int n) : n(n), adj(n), first(n), level(n), weight_sum(n) {}

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

    long long dist(int u, int v) const{
        return weight_sum[u] + weight_sum[v] - 2 * weight_sum[lca(u, v)];
    }

    int lca(int u, int v) const{
        int l = first[u], r = first[v];
        if (l > r) swap(l, r);
        return tour[rmq.query(l, r)];
    }
};
/// END COPY linear_lca

struct VirtualTree : LinearLCA{
    using LinearLCA::LinearLCA;

    struct Node{ int v, parent; long long len; };

    /// The parent of each vertex is the lca of it and its predecessor in DFS order, so no stack is needed
    vector<Node> compress(vector<int> nodes) const{
        auto by_tour = [&](int u, int v){ return first[u] < first[v]; };

        sort(nodes.begin(), nodes.end(), by_tour);
        int k = nodes.size();
        for (int i = 0; i + 1 < k; i++) nodes.push_back(lca(nodes[i], nodes[i + 1]));
        sort(nodes.begin(), nodes.end(), by_tour);
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.empty()) return {};

        vector<Node> res = {{nodes[0], -1, 0}};
        for (int i = 1; i < (int)nodes.size(); i++){
            int p = lca(nodes[i - 1], nodes[i]);
            int p_idx = lower_bound(nodes.begin(), nodes.end(), p, by_tour) - nodes.begin();
            res.push_back({nodes[i], p_idx, weight_sum[nodes[i]] - weight_sum[p]});
        }

        return res;
    }
};

int main(){
    /***
     *             0
     *        3 /     \ 1
     *         1       2
     *    2 /    \ 7    \ 4
     *     3      4      5
     *  5 /
     *   6
    ***/
    vector<array<int, 3>> edges = {{0, 1, 3}, {0, 2, 1}, {1, 3, 2}, {1, 4, 7}, {2, 5, 4}, {3, 6, 5}};
    VirtualTree tree(7);
    for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
    tree.build(0);

    using Result = vector<array<long long, 3>>;
    auto fields = [](const vector<VirtualTree::Node>& vt){
        Result res;
        for (auto [v, parent, len] : vt) res.push_back({v, parent, len});
        return res;
    };

    /// 6 and 4 meet at 1, which meets 5 at 0; compressed edge 6 - 1 stands for 6 3 1, 5 + 2
    assert((fields(tree.compress({6, 4, 5})) == Result{{0, -1, 0}, {1, 0, 3}, {6, 1, 7}, {4, 1, 7}, {5, 0, 5}}));
    assert((fields(tree.compress({6, 3})) == Result{{3, -1, 0}, {6, 0, 5}}));
    assert(tree.compress({}).empty());
    assert(tree.lca(6, 4) == 1);
    assert(tree.depth(6) == 3);
    return 0;
}
