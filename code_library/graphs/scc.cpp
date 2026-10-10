/***
 *
 * Strongly Connected Components (Kosaraju)
 * Groups the vertices of a directed graph into strongly connected components
 *
 * Complexity: O(n + m)
 *
 * SCC g(n); g.add_edge(u, v); g.run();
 * g.count: number of components, g.comp[v]: component of v in [0, count)
 * Components come out in topological order of the condensation: every edge u -> v has comp[u] <= comp[v]
 * g.condensation(): adjacency lists of the component DAG, without duplicate edges or self loops
 *
 * Both passes are iterative, so long paths cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SCC{
    int n, count = 0;
    vector<vector<int>> adj, rev;
    vector<int> comp;

    SCC(int n) : n(n), adj(n), rev(n), comp(n, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        rev[v].push_back(u);
    }

    void run(){
        vector<int> order, it(n, 0), stack;
        vector<char> visited(n, 0);
        for (int s = 0; s < n; s++){
            if (visited[s]) continue;
            stack.push_back(s);
            visited[s] = 1;
            while (!stack.empty()){
                int u = stack.back();
                if (it[u] < (int)adj[u].size()){
                    int v = adj[u][it[u]++];
                    if (!visited[v]) visited[v] = 1, stack.push_back(v);
                }
                else{
                    order.push_back(u);
                    stack.pop_back();
                }
            }
        }

        count = 0;
        fill(comp.begin(), comp.end(), -1);
        for (int i = n - 1; i >= 0; i--){
            int s = order[i];
            if (comp[s] != -1) continue;
            stack.push_back(s);
            comp[s] = count;
            while (!stack.empty()){
                int u = stack.back();
                stack.pop_back();
                for (int v : rev[u]){
                    if (comp[v] == -1) comp[v] = count, stack.push_back(v);
                }
            }
            count++;
        }
    }

    vector<vector<int>> condensation() const{
        vector<vector<int>> dag(count);
        for (int u = 0; u < n; u++){
            for (int v : adj[u]){
                if (comp[u] != comp[v]) dag[comp[u]].push_back(comp[v]);
            }
        }
        for (auto& out : dag){
            sort(out.begin(), out.end());
            out.erase(unique(out.begin(), out.end()), out.end());
        }
        return dag;
    }
};

int main(){
    /***
     * 0 -> 1 -> 2 -> 0 is a cycle, 2 -> 3, 3 -> 4 -> 3 is a cycle, 5 is alone with a self loop
    ***/
    SCC g(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 3}, {5, 5}}) g.add_edge(u, v);
    g.run();
    assert(g.count == 3);
    assert(g.comp[0] == g.comp[1] && g.comp[1] == g.comp[2]);
    assert(g.comp[3] == g.comp[4]);
    assert(g.comp[0] < g.comp[3]);  /// topological order, edge 2 -> 3

    auto dag = g.condensation();
    assert((dag[g.comp[0]] == vector<int>{g.comp[3]}));
    assert(dag[g.comp[5]].empty());  /// the self loop is dropped
    return 0;
}
