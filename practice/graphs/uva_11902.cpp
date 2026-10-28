/***
 *
 * Dominator Tree (Lengauer-Tarjan)
 * Immediate dominator of every vertex of a directed graph, for a fixed root
 *
 * Complexity: O((n + m) log n), path compression without union by rank
 *
 * d dominates v if every path from root to v passes through d, idom[v] is the closest strict dominator
 * The idom array is the dominator tree as a parent array: d dominates v iff d is an ancestor of v in it
 *
 * DominatorTree g(n); g.add_edge(u, v); auto idom = g.run(root);
 * idom[root] = root, idom[v] = -1 if v is unreachable from root
 * Vertices are 0-indexed, parallel edges and self loops allowed, run can be called again with another root
 *
 * The DFS and the path compression are iterative, so long paths cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DominatorTree{
    int n;
    vector<vector<int>> adj;

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

    vector<int> order, num, parent, sdom, label, anc, path;

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
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);
        DominatorTree g(n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                int x;
                scanf("%d", &x);
                if (x) g.add_edge(i, j);
            }
        }

        vector<int> idom = g.run(0);
        vector<string> grid(n, string(n, 'N'));
        for (int v = 0; v < n; v++){
            if (idom[v] == -1) continue;
            for (int d = v; ; d = idom[d]){
                grid[d][v] = 'Y';
                if (d == 0) break;
            }
        }

        string border = "+" + string(2 * n - 1, '-') + "+";
        printf("Case %d:\n%s\n", tc, border.c_str());
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) printf("|%c", grid[i][j]);
            printf("|\n%s\n", border.c_str());
        }
    }

    return 0;
}
