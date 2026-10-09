/***
 *
 * Finds articulation points in simple graphs
 *
 * An articulation point or cut vertex is defined as a vertex which, when removed along with associated edges, makes the graph disconnected
 * More precisely, increases the number of connected components in the graph
 *
 * Nodes are numbered from 0 to N-1
 * Complexity: O(N + M)
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Graph{
    int n, disc_t;
    vector <vector<int>> adj;
    vector <int> discover, low;
    vector <bool> is_cut;

    Graph(int n) : n(n), adj(n) {}

    void dfs(int u, int p){
        discover[u] = low[u] = ++disc_t;

        int children = 0;
        for (auto v: adj[u]){
            if (!discover[v]){
                children++;
                dfs(v, u);
                low[u] = min(low[u], low[v]);

                if ((low[v] >= discover[u]) && !(p == -1 && children < 2)) {
                    is_cut[u] = true;
                }
            }
            else if (v != p) low[u] = min(low[u], discover[v]);
        }
    }

    /// adds undirected edge from u to v
    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector <int> get_cuts(){
        discover.assign(n, 0), low.assign(n, 0), is_cut.assign(n, false);

        vector <int> cuts;
        for (int i = 0; i < n; i++){
            if (!discover[i]){
                disc_t = 0;
                dfs(i, -1);
            }
            if (is_cut[i]) cuts.push_back(i);
        }

        return cuts;
    }
};

int main(){
    auto graph = Graph(10);

    graph.add_edge(0, 1);
    graph.add_edge(1, 2);
    graph.add_edge(2, 0);
    graph.add_edge(5, 6);
    graph.add_edge(6, 7);
    graph.add_edge(7, 6);
    graph.add_edge(2, 3);
    graph.add_edge(3, 4);
    graph.add_edge(2, 4);
    graph.add_edge(4, 5);
    graph.add_edge(8, 9);

    auto cuts = graph.get_cuts();
    assert((int)cuts.size() == 4);
    assert(cuts == vector<int>({2, 4, 5, 6}));

    return 0;
}
