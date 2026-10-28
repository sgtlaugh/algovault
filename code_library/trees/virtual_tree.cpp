/***
 *
 * Virtual Tree (Auxiliary Tree)
 * Compresses a rooted tree to a vertex subset plus the pairwise LCAs needed to keep its shape, 0-based nodes
 *
 * Complexity: O(n log n) to build and memory, O(k log k) per compress of k vertices
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

/// BEGIN COPY fast_lca from code_library/trees/lca.cpp
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
/// END COPY fast_lca

struct VirtualTree : FastLCA{
    using FastLCA::FastLCA;

    struct Node{ int v, parent; long long len; };

    /// The parent of each vertex is the lca of it and its predecessor in DFS order, so no stack is needed
    vector<Node> compress(vector<int> nodes) const{
        auto by_tin = [&](int u, int v){ return tin[u] < tin[v]; };

        sort(nodes.begin(), nodes.end(), by_tin);
        int k = nodes.size();
        for (int i = 0; i + 1 < k; i++) nodes.push_back(lca(nodes[i], nodes[i + 1]));
        sort(nodes.begin(), nodes.end(), by_tin);
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        if (nodes.empty()) return {};

        vector<Node> res = {{nodes[0], -1, 0}};
        for (int i = 1; i < (int)nodes.size(); i++){
            int p = lca(nodes[i - 1], nodes[i]);
            int p_idx = lower_bound(nodes.begin(), nodes.end(), p, by_tin) - nodes.begin();
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
    assert((fields(tree.compress({6, 4, 5})) == Result{{0, -1, 0}, {5, 0, 5}, {1, 0, 3}, {4, 2, 7}, {6, 2, 7}}));
    assert((fields(tree.compress({6, 3})) == Result{{3, -1, 0}, {6, 0, 5}}));
    assert(tree.compress({}).empty());
    assert(tree.lca(6, 4) == 1);
    assert(tree.depth(6) == 3);
    return 0;
}
