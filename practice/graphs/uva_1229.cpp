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
    int n;
    static char line[1 << 12];
    while (scanf("%d", &n) == 1 && n){
        fgets(line, sizeof(line), stdin);
        map<string, int> id;
        auto get = [&](const string& s){ return id.emplace(s, id.size()).first->second; };
        vector<pair<int, int>> edges;
        for (int i = 0; i < n; i++){
            fgets(line, sizeof(line), stdin);
            istringstream in(line);
            string word, def;
            in >> word;
            int u = get(word);
            while (in >> def) edges.push_back({u, get(def)});
        }

        SCC g(n);
        vector<int> size(n, 0);
        for (auto edge : edges) g.add_edge(edge.first, edge.second);
        g.run();
        for (int v = 0; v < n; v++) size[g.comp[v]]++;

        /// A word on a cycle can never be learned, so the sub-dictionary must hold every cycle and, being a dictionary
        /// itself, everything their definitions reach
        vector<char> chosen(n, 0);
        vector<int> stack;
        for (int v = 0; v < n; v++){
            if (size[g.comp[v]] > 1) chosen[v] = 1, stack.push_back(v);
        }
        while (!stack.empty()){
            int u = stack.back();
            stack.pop_back();
            for (int v : g.adj[u]){
                if (!chosen[v]) chosen[v] = 1, stack.push_back(v);
            }
        }

        vector<string> res;
        for (auto& entry : id){
            if (chosen[entry.second]) res.push_back(entry.first);
        }
        printf("%d\n", (int)res.size());
        for (int i = 0; i < (int)res.size(); i++) printf("%s%s", i ? " " : "", res[i].c_str());
        puts("");
    }
    return 0;
}
