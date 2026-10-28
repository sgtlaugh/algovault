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
    int t, n;

    scanf("%d", &t);
    while (t--){
        scanf("%d", &n);
        vector<long long> x(n), y(n), d(n);
        for (int i = 0; i < n; i++) scanf("%lld %lld %lld", &x[i], &y[i], &d[i]);

        SCC g(n);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (i != j && 2 * abs(x[i] - x[j]) <= d[i] && 2 * abs(y[i] - y[j]) <= d[i]) g.add_edge(i, j);
            }
        }
        g.run();

        vector<bool> reached(g.count, false);
        for (auto& out : g.condensation()){
            for (int c : out) reached[c] = true;
        }
        printf("%d\n", (int)count(reached.begin(), reached.end(), false));
    }

    return 0;
}
