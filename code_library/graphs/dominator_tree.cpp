/***
 *
 * Dominator Tree (Lengauer-Tarjan)
 * Immediate dominator of every vertex of a directed graph, for a fixed root
 *
 * Complexity: O((n + m) log n)
 *
 * DominatorTree g(n): vertices are 0 to n - 1, parallel edges and self loops allowed
 * g.add_edge(u, v): directed edge u -> v
 * g.run(root): idom, idom[v] the closest strict dominator of v, idom[root] = root, -1 if v is unreachable from root
 *
 * Notes:
 *   - d dominates v if every path from root to v passes through d
 *   - idom is the dominator tree as a parent array: d dominates v iff d is an ancestor of v in it
 *   - run can be called again with another root
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DominatorTree{
    int n;
    vector<vector<int>> adj;
    vector<int> order, num, parent, sdom, label, anc, path;

    DominatorTree(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    vector<int> run(int root){
        dfs(root);
        int k = order.size();
        vector<vector<int>> pred(k), bucket(k);
        for (int i = 0; i < k; i++){
            for (int v : adj[order[i]]) pred[num[v]].push_back(i);
        }

        vector<int> dom(k, 0);
        sdom.resize(k), label.resize(k), anc.resize(k);
        for (int i = 0; i < k; i++) sdom[i] = label[i] = anc[i] = i;
        for (int i = k - 1; i >= 0; i--){
            for (int p : pred[i]) sdom[i] = min(sdom[i], sdom[eval(p)]);
            if (i > 0) bucket[sdom[i]].push_back(i);
            for (int w : bucket[i]){
                int v = eval(w);
                dom[w] = sdom[v] == sdom[w] ? sdom[w] : v;
            }
            if (i > 0) anc[i] = parent[i];
        }

        vector<int> idom(n, -1);
        idom[root] = root;
        for (int i = 1; i < k; i++){
            if (dom[i] != sdom[i]) dom[i] = dom[dom[i]];
            idom[order[i]] = order[dom[i]];
        }
        return idom;
    }

    /// Preorder numbering, parent[i] is the DFS tree parent of the i-th visited vertex, num[v] = -1 if unreachable
    void dfs(int root){
        order.clear(), parent.clear();
        num.assign(n, -1);
        vector<int> it(n, 0), stack = {root};
        num[root] = 0, order.push_back(root), parent.push_back(0);
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == (int)adj[u].size()){
                stack.pop_back();
                continue;
            }

            int v = adj[u][it[u]++];
            if (num[v] != -1) continue;
            num[v] = order.size();
            order.push_back(v), parent.push_back(num[u]);
            stack.push_back(v);
        }
    }

    /// Vertex of minimum sdom on the forest path from v up to its root, the root itself excluded
    int eval(int v){
        if (anc[v] == v) return v;

        path.clear();
        for (int u = v; anc[anc[u]] != anc[u]; u = anc[u]) path.push_back(u);
        for (int i = (int)path.size() - 1; i >= 0; i--){
            int w = path[i], a = anc[w];
            if (sdom[label[a]] < sdom[label[w]]) label[w] = label[a];
            anc[w] = anc[a];
        }
        return label[v];
    }
};

int main(){
    /// Figure 1 of Lengauer and Tarjan, "A Fast Algorithm for Finding Dominators in a Flowgraph", R to L are 0 to 12
    vector<vector<int>> out = {
        {1, 2, 3}, {4}, {1, 4, 5}, {6, 7}, {12}, {8}, {9}, {9, 10}, {5, 11}, {11}, {9}, {9, 0}, {8}
    };
    DominatorTree paper(13);
    for (int u = 0; u < 13; u++){
        for (int v : out[u]) paper.add_edge(u, v);
    }
    assert((paper.run(0) == vector<int>{0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 7, 0, 4}));

    /// 0 reaches 3 through 1 or 2, so only 0 dominates 3, and 5 never gets reached from 0
    vector<pair<int, int>> edges = {{0, 1}, {0, 2}, {1, 3}, {2, 3}, {3, 4}, {4, 4}, {5, 0}};
    DominatorTree diamond(6);
    for (auto [u, v] : edges) diamond.add_edge(u, v);
    assert((diamond.run(0) == vector<int>{0, 0, 0, 0, 3, -1}));
    assert((diamond.run(5) == vector<int>{5, 0, 0, 0, 3, 5}));  /// the same graph from another root

    return 0;
}
