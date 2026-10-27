/***
 *
 * DSU on Tree (small to large, sack)
 * Answers a query about every subtree by keeping the heavy child's data and re-adding only light subtrees
 *
 * Complexity: O(n log n) calls to add / remove
 *
 * DsuOnTree tree(n): an empty tree on nodes 0..n-1
 * tree.add_edge(u, v): adds the edge u - v
 * tree.run(root, add, remove, answer): calls answer(v) for every v while the structure holds exactly v's subtree
 * tree.tin[v]: the position of v in order
 * tree.tout[v]: the last position of v's subtree, so the subtree of v is order[tin[v] .. tout[v]]
 * tree.order[i]: the node at position i of a preorder that visits each node's heavy child first
 * tree.parent[v]: the parent of v, -1 for the root
 * tree.size[v]: the number of nodes in the subtree of v
 * tree.heavy[v]: a child of v with the largest subtree, -1 for a leaf
 *
 * Notes:
 *   - The fields describe the tree rooted at the root of the last run call
 *   - add(v) / remove(v) put node v into / take it out of your data structure, empty again when run returns
 *   - The graph must be a tree, run asserts that every node is reached exactly once
 *   - Iterative, so deep trees cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DsuOnTree{
    int n;
    vector<vector<int>> adj;
    vector<int> parent, heavy, size, tin, tout, order;

    DsuOnTree(int n) : n(n), adj(n), parent(n, -1), heavy(n, -1), size(n, 1), tin(n), tout(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    template <typename Add, typename Remove, typename Answer>
    void run(int root, Add add, Remove remove, Answer answer){
        prepare(root);

        /// Reverse preorder reaches v right after its heavy child, whose subtree the structure still holds
        for (int i = n - 1; i >= 0; i--){
            int v = order[i];
            add(v);
            for (int j = heavy[v] == -1 ? i + 1 : tout[heavy[v]] + 1; j <= tout[v]; j++) add(order[j]);
            answer(v);
            if (v == root || heavy[parent[v]] != v){
                for (int j = i; j <= tout[v]; j++) remove(order[j]);
            }
        }
    }

    /// Preorder with the heavy child first, so every subtree is one contiguous range
    void prepare(int root){
        size.assign(n, 1), heavy.assign(n, -1), parent.assign(n, -1);  /// run may be called again, from any root
        vector<int> bfs = {root};
        parent[root] = -1;
        for (int i = 0; i < (int)bfs.size(); i++){
            for (int v : adj[bfs[i]]){
                if (v == parent[bfs[i]]) continue;
                assert(v != root && parent[v] == -1);  /// a revisit means a cycle, which would grow bfs forever
                parent[v] = bfs[i], bfs.push_back(v);
            }
        }
        assert((int)bfs.size() == n);

        for (int i = n - 1; i > 0; i--){
            int v = bfs[i], p = parent[v];
            size[p] += size[v];
            if (heavy[p] == -1 || size[v] > size[heavy[p]]) heavy[p] = v;
        }

        order.clear();
        vector<int> stack = {root};
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            tin[u] = order.size(), order.push_back(u);
            for (int v : adj[u]){
                if (v != parent[u] && v != heavy[u]) stack.push_back(v);
            }
            if (heavy[u] != -1) stack.push_back(heavy[u]);
        }

        for (int i = n - 1; i >= 0; i--) tout[order[i]] = tin[order[i]] + size[order[i]] - 1;
    }
};

int main(){
    /***
     *          0 (color 1)
     *        /   \
     *       1(2)  2(1)
     *      / \     \
     *    3(2) 4(3)  5(2)
     * Distinct colors per subtree: 0 -> 3, 1 -> 2, 2 -> 2, leaves -> 1
    ***/
    DsuOnTree tree(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}}) tree.add_edge(u, v);
    vector<int> color = {1, 2, 1, 2, 3, 2}, freq(4, 0), distinct(6, -1);
    int count = 0;
    tree.run(0,
        [&](int v){ count += freq[color[v]]++ == 0; },
        [&](int v){ count -= --freq[color[v]] == 0; },
        [&](int v){ distinct[v] = count; });

    assert((distinct == vector<int>{3, 2, 2, 1, 1, 1}));
    assert(count == 0);

    assert(tree.tout[1] - tree.tin[1] == 2);    /// order[tin[1] .. tout[1]] holds 1, 3, 4
    assert(tree.order[tree.tin[4]] == 4);
    return 0;
}
