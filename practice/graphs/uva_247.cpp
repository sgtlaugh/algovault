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
    int n, m;

    for (int tc = 1; scanf("%d %d", &n, &m) == 2 && (n || m); tc++){
        map<string, int> id;
        vector<string> names;
        vector<pair<int, int>> calls;
        auto get = [&](const string& s){
            auto it = id.find(s);
            if (it != id.end()) return it->second;
            names.push_back(s);
            return id[s] = names.size() - 1;
        };

        char a[32], b[32];
        for (int i = 0; i < m; i++){
            scanf("%s %s", a, b);
            int u = get(a);
            calls.push_back({u, get(b)});
        }

        SCC g(names.size());
        for (auto call : calls) g.add_edge(call.first, call.second);
        g.run();

        vector<vector<int>> circles(g.count);
        vector<int> order;
        for (int v = 0; v < (int)names.size(); v++){
            if (circles[g.comp[v]].empty()) order.push_back(g.comp[v]);
            circles[g.comp[v]].push_back(v);
        }

        if (tc > 1) puts("");
        printf("Calling circles for data set %d:\n", tc);
        for (int c : order){
            for (int i = 0; i < (int)circles[c].size(); i++){
                printf("%s%s", i ? ", " : "", names[circles[c][i]].c_str());
            }
            puts("");
        }
    }

    return 0;
}
