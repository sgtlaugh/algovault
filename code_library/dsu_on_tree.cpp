/***
 *
 * DSU on Tree (small to large, sack)
 * Answers a query about every subtree by keeping the heavy child's data and re-adding only light subtrees
 *
 * Complexity: O(n log n) calls to add / remove
 *
 * DsuOnTree tree(n); tree.add_edge(u, v); tree.run(root, add, remove, answer);
 *   add(v) / remove(v) put node v into / take it out of your data structure
 *   answer(v) is called when the structure holds exactly the nodes of v's subtree
 *   the structure is empty again when run returns
 * Subtree positions (Euler tour), filled by run:
 *   tree.tin[v], tree.tout[v]: the subtree of v is order[tin[v] .. tout[v]], tree.order[i]: node at position i
 *
 * Iterative, so deep trees cannot overflow the stack
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

    /// Preorder with the heavy child first, so every subtree is one contiguous range
    void prepare(int root){
        vector<int> bfs = {root};
        parent[root] = -1;
        for (int i = 0; i < (int)bfs.size(); i++){
            for (int v : adj[bfs[i]]){
                if (v != parent[bfs[i]]) parent[v] = bfs[i], bfs.push_back(v);
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

    template <typename Add, typename Remove, typename Answer>
    void run(int root, Add add, Remove remove, Answer answer){
        prepare(root);

        /// Frame (v, keep, stage): stage 0 = solve light children, 1 = solve heavy child, 2 = add the rest and answer
        vector<tuple<int, bool, int>> stack = {{root, true, 0}};
        while (!stack.empty()){
            auto [v, keep, stage] = stack.back();  /// a copy: pushing below may reallocate the stack
            if (stage == 0){
                get<2>(stack.back()) = 1;
                for (int c : adj[v]){
                    if (c != parent[v] && c != heavy[v]) stack.push_back({c, false, 0});
                }
            }
            else if (stage == 1){
                get<2>(stack.back()) = 2;
                if (heavy[v] != -1) stack.push_back({heavy[v], true, 0});
            }
            else{
                stack.pop_back();
                add(v);
                for (int c : adj[v]){
                    if (c == parent[v] || c == heavy[v]) continue;
                    for (int i = tin[c]; i <= tout[c]; i++) add(order[i]);
                }
                answer(v);
                if (!keep){
                    for (int i = tin[v]; i <= tout[v]; i++) remove(order[i]);
                }
            }
        }
        for (int i = 0; i < n; i++) remove(order[i]);
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
    assert(tree.tout[1] - tree.tin[1] == 2 && tree.tout[0] - tree.tin[0] == 5 && tree.tout[2] - tree.tin[2] == 1);
    assert(tree.order[tree.tin[4]] == 4);

    DsuOnTree single(1);
    int calls = 0;
    single.run(0, [&](int){}, [&](int){}, [&](int){ calls++; });
    assert(calls == 1);
    return 0;
}
