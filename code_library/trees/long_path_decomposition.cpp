/***
 *
 * Long-Path Decomposition
 * Splits a rooted tree into vertex-disjoint paths that always continue into the child with the deepest subtree
 *
 * Complexity: O(n log n) build and memory, O(1) per kth_ancestor, O(n) per depth_counts pass plus the cost of visit
 *
 * LongPathDecomposition tree(n); tree.add_edge(u, v); tree.build(root);   n >= 1, 0-based nodes, unweighted edges
 *   - kth_ancestor(v, k): the k-th ancestor of v, v itself for k = 0, -1 if k < 0 or k > depth(v)
 *     Ladder algorithm: every path is extended upwards by its own length, one binary lifting jump by the highest
 *     bit of k lands on a node whose path's ladder reaches the answer, so the rest is an array lookup
 *   - depth_counts(visit): calls visit(v, cnt, len) once per node, children before parents, where
 *     cnt[d] for 0 <= d < len = height(v) + 1 is the number of nodes in v's subtree exactly d edges below v
 *     cnt points into a shared buffer that later merges overwrite, so read it during the call only
 *     A node inherits its long child's array in place and merges only the other children, which is O(n) in total
 *   - depth(v): edges from the root to v, height(v): edges on the longest downward path from v
 *   The graph must be a tree, build asserts n - 1 edges and that every node is reached exactly once
 *   Traversals are iterative so deep trees are safe, build may be called again with another root
 *
 * Offline "how many nodes are d edges below v": group the queries by v, then inside visit(v, cnt, len)
 *   each query (v, d) is d < len ? cnt[d] : 0
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct LongPathDecomposition{
    int n, lg, added_edges = 0;
    vector<vector<int>> adj, up;
    vector<int> order, parent, level, heights, long_child, pos, ladder, ladder_index;

    LongPathDecomposition(int n) : n(n), adj(n), level(n), pos(n), ladder_index(n){
        lg = 1;
        while ((1 << lg) < n) lg++;
    }

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
        added_edges++;
    }

    void build(int root = 0){
        assert(added_edges == n - 1);
        parent.assign(n, -1);
        order = {root};
        parent[root] = root, level[root] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int v : adj[u]){
                if (v == parent[u]) continue;
                assert(parent[v] == -1);
                parent[v] = u, level[v] = level[u] + 1;
                order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        heights.assign(n, 0), long_child.assign(n, -1);
        for (int i = n - 1; i > 0; i--){
            int v = order[i], p = parent[v];
            if (heights[v] + 1 > heights[p]) heights[p] = heights[v] + 1, long_child[p] = v;
        }

        /// ladder of a path: up to len ancestors of its head, then the path itself, at most 2n entries in total
        ladder.clear(), ladder.reserve(2 * n);
        int next_pos = 0;
        for (int head : order){
            if (head != root && long_child[parent[head]] == head) continue;
            int len = heights[head] + 1, above = min(len, level[head]), start = ladder.size();
            ladder.resize(start + above + len);
            for (int j = 1, a = head; j <= above; j++){
                a = parent[a];
                ladder[start + above - j] = a;
            }
            for (int i = 0, v = head; i < len; i++, v = long_child[v]){
                pos[v] = next_pos++, ladder_index[v] = start + above + i;
                ladder[ladder_index[v]] = v;
            }
        }

        up.assign(lg, parent);
        for (int k = 1; k < lg; k++){
            for (int v = 0; v < n; v++) up[k][v] = up[k - 1][up[k - 1][v]];
        }
    }

    int depth(int v) const{
        return level[v];
    }

    /// Path nodes have consecutive positions, so a node's array is its long child's array shifted down by one slot
    template <typename Visit>
    void depth_counts(Visit visit) const{
        vector<int> cnt(n, 0);
        for (int i = n - 1; i >= 0; i--){
            int v = order[i];
            cnt[pos[v]] = 1;
            for (int c : adj[v]){
                if (c == parent[v] || c == long_child[v]) continue;
                for (int d = 0; d <= heights[c]; d++) cnt[pos[v] + 1 + d] += cnt[pos[c] + d];
            }
            visit(v, (const int*)cnt.data() + pos[v], heights[v] + 1);
        }
    }

    int height(int v) const{
        return heights[v];
    }

    /// The jump lands on u with height(u) >= 2^h > the remaining distance, and u's path is at least that tall, so its ladder covers it
    int kth_ancestor(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        if (k == 0) return v;

        int h = __lg(k), u = up[h][v];
        return ladder[ladder_index[u] - (k - (1 << h))];
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
    LongPathDecomposition tree(7);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {3, 6}}) tree.add_edge(u, v);
    tree.build(0);

    vector<vector<int>> counts(7);
    auto collect = [&](int v, const int* cnt, int len){ counts[v].assign(cnt, cnt + len); };
    tree.depth_counts(collect);
    assert((counts == vector<vector<int>>{{1, 2, 3, 1}, {1, 2, 1}, {1, 1}, {1, 1}, {1}, {1}, {1}}));

    assert(tree.depth(6) == 3 && tree.depth(0) == 0 && tree.height(0) == 3 && tree.height(1) == 2 && tree.height(4) == 0);
    assert(tree.kth_ancestor(6, 0) == 6);
    assert(tree.kth_ancestor(6, 1) == 3);
    assert(tree.kth_ancestor(6, 2) == 1);
    assert(tree.kth_ancestor(6, 3) == 0);
    assert(tree.kth_ancestor(6, 4) == -1);
    assert(tree.kth_ancestor(5, 2) == 0);
    assert(tree.kth_ancestor(4, 1) == 1);
    assert(tree.kth_ancestor(2, -1) == -1);

    tree.build(4);
    tree.depth_counts(collect);
    assert((counts[4] == vector<int>{1, 1, 2, 2, 1}));
    assert((counts[1] == vector<int>{1, 2, 2, 1}));
    assert((counts[0] == vector<int>{1, 1, 1}));
    assert(tree.kth_ancestor(5, 4) == 4);
    assert(tree.kth_ancestor(5, 3) == 1);
    assert(tree.kth_ancestor(6, 2) == 1);
    assert(tree.depth(5) == 4 && tree.height(4) == 4 && tree.height(3) == 1);

    LongPathDecomposition path(10);
    for (int i = 0; i + 1 < 10; i++) path.add_edge(i, i + 1);
    path.build(0);
    for (int k = 0; k <= 9; k++) assert(path.kth_ancestor(9, k) == 9 - k);
    assert(path.kth_ancestor(9, 10) == -1);

    LongPathDecomposition single(1);
    single.build(0);
    int calls = 0;
    single.depth_counts([&](int v, const int* cnt, int len){ calls++, assert(v == 0 && len == 1 && cnt[0] == 1); });
    assert(calls == 1);
    assert(single.kth_ancestor(0, 0) == 0 && single.kth_ancestor(0, 1) == -1);

    return 0;
}
