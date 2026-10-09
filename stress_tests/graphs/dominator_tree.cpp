#include "../common.h"

#define main library_main
#include "../../code_library/graphs/dominator_tree.cpp"
#undef main

/// Reachable vertices from root, skipping the vertex removed (-1 removes none)
vector<char> reach(const vector<vector<int>>& adj, int root, int removed){
    vector<char> seen(adj.size(), 0);
    if (root == removed) return seen;
    vector<int> queue = {root};
    seen[root] = 1;
    for (int i = 0; i < (int)queue.size(); i++){
        for (int v : adj[queue[i]]){
            if (v != removed && !seen[v]) seen[v] = 1, queue.push_back(v);
        }
    }
    return seen;
}

/// d strictly dominates v iff removing d cuts v off, idom is the strict dominator with the most strict dominators of its own
vector<int> brute_idom(int n, const vector<pair<int, int>>& edges, int root){
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges) adj[u].push_back(v);

    vector<char> base = reach(adj, root, -1);
    vector<vector<char>> dom(n, vector<char>(n, 0));
    vector<int> count(n, 0);
    for (int d = 0; d < n; d++){
        vector<char> cut = reach(adj, root, d);
        for (int v = 0; v < n; v++){
            if (v != d && base[v] && !cut[v]) dom[d][v] = 1, count[v]++;
        }
    }

    vector<int> idom(n, -1);
    idom[root] = root;
    for (int v = 0; v < n; v++){
        if (!base[v] || v == root) continue;
        idom[v] = root;
        for (int d = 0; d < n; d++){
            if (dom[d][v] && count[d] > count[idom[v]]) idom[v] = d;
        }
    }
    return idom;
}

/// One instance answers every root, so state left over from an earlier run is under the oracle too
void check(int n, const vector<pair<int, int>>& edges, const vector<int>& roots){
    DominatorTree g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);
    for (int root : roots) assert(g.run(root) == brute_idom(n, edges, root));
}

vector<int> random_roots(int n, int count){
    vector<int> roots;
    for (int i = 0; i < count; i++) roots.push_back(stress::rand_int(0, n - 1));
    return roots;
}

int main(){
    /// Every graph on 3 vertices with self loops, every graph on 4 vertices without, every root
    for (int n = 3; n <= 4; n++){
        vector<int> roots(n);
        iota(roots.begin(), roots.end(), 0);
        vector<pair<int, int>> all;
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++){
                if (n == 3 || u != v) all.push_back({u, v});
            }
        }
        for (int mask = 0; mask < (1 << all.size()); mask++){
            vector<pair<int, int>> edges;
            for (int i = 0; i < (int)all.size(); i++){
                if (mask >> i & 1) edges.push_back(all[i]);
            }
            check(n, edges, roots);
        }
    }

    /// Sparse graphs give long dominator chains, dense ones collapse onto the root
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, it < 3000 ? 12 : 120);
        int m = stress::rand_int(0, stress::rand_int(0, 3) ? 2 * n : n * n / 4);
        vector<pair<int, int>> edges;
        for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        check(n, edges, random_roots(n, 3));
    }

    /// Random forward tree plus back and cross edges, so DFS order and dominators disagree often
    for (long long it = 0; it < stress::scaled(10); it++){
        int n = stress::rand_int(500, 1500);
        vector<pair<int, int>> edges;
        for (int v = 1; v < n; v++) edges.push_back({(int)stress::rand_int(max(0, v - 5), v - 1), v});
        for (int i = 0; i < n; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        check(n, edges, {0, (int)stress::rand_int(1, n - 1), 0});
    }

    /// A path deep enough to overflow a recursive DFS or a recursive eval on an 8 MB stack, plus a back edge to every vertex from the end
    int n = 300000;
    DominatorTree path(n), ladder(n);
    for (int i = 0; i + 1 < n; i++) path.add_edge(i, i + 1), ladder.add_edge(i, i + 1);
    for (int i = 1; i + 1 < n; i++) ladder.add_edge(n - 1, i), ladder.add_edge(0, i);

    vector<int> chain = path.run(0), flat = ladder.run(0);
    assert(chain[0] == 0 && flat[0] == 0);
    for (int i = 1; i < n; i++) assert(chain[i] == i - 1 && flat[i] == (i == n - 1 ? n - 2 : 0));
    return 0;
}
