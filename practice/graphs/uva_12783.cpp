/***
 *
 * Bridges, 2-edge-connected components and the bridge tree of an undirected graph, parallel edges and self loops allowed
 * A doubled edge is never a bridge: the search skips only the edge it arrived by, not every edge to the parent
 *
 * A bridge is such an edge which when removed makes the graph disconnected
 * Or more precisely, increases the number of connected components
 *
 * BridgeGraph g(n); int id = g.add_edge(u, v); nodes are numbered from 0 to n - 1
 * g.get_bridges(): every bridge with its endpoints, its edge id and the node counts on both sides, so weights or labels can be kept outside
 * g.get_components(): groups of the 2-edge-connected components, the components left after removing every bridge
 *     also sets g.comp[v] in [0, g.count), numbered by smallest vertex, each group ascending
 * g.get_bridge_tree(): adjacency lists over the g.count components, one tree edge per bridge, see
 *     https://www.quora.com/q/threadsiiithyderabad/The-Bridge-Tree-of-a-graph
 *     a forest, one tree per connected component
 *
 * Both searches are iterative with heap-allocated stacks, so a 2e5-node path needs no extra call stack
 *
 * Complexity: O(n + m)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Bridge{
    int u, v; /// bridge edge from node u to v
    int cnt_u, cnt_v; /// number of nodes in the connected component of u and v if bridge edge is disconnected
    int id; /// index of the edge, in the order add_edge was called

    Bridge() {}
    Bridge(int u, int v, int cnt_u, int cnt_v, int id = -1) : u(u), v(v), cnt_u(cnt_u), cnt_v(cnt_v), id(id) {}
};

struct BridgeGraph{
    int n, m = 0, count = 0;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge index)
    vector<int> comp;

    BridgeGraph(int n): n(n), adj(n), comp(n, -1) {}

    int add_edge(int u, int v){
        adj[u].push_back({v, m});
        adj[v].push_back({u, m});
        return m++;
    }

    /// Bridges are listed in the order the DFS finishes their endpoint v
    vector<Bridge> get_bridges(){
        vector<Bridge> bridges;
        vector<int> discover(n, 0), low(n), parent_edge(n), edge_pos(n, 0), dfs_stack;

        for (int i = 0; i < n; i++){
            if (discover[i]) continue;

            int first = bridges.size(), dt = 0;
            discover[i] = low[i] = ++dt, parent_edge[i] = -1;
            dfs_stack.push_back(i);

            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                if (edge_pos[u] < (int)adj[u].size()){
                    int v, id;
                    tie(v, id) = adj[u][edge_pos[u]++];
                    if (!discover[v]){
                        discover[v] = low[v] = ++dt, parent_edge[v] = id;
                        dfs_stack.push_back(v);
                    }
                    else if (id != parent_edge[u]) low[u] = min(low[u], discover[v]);
                    continue;
                }

                dfs_stack.pop_back();
                if (dfs_stack.empty()) break;

                int p = dfs_stack.back();
                low[p] = min(low[p], low[u]);
                if (low[u] > discover[p]) bridges.emplace_back(p, u, 0, dt - discover[u] + 1, parent_edge[u]);
            }

            /// The component size is known only once its DFS ends
            for (int j = first; j < (int)bridges.size(); j++) bridges[j].cnt_u = dt - bridges[j].cnt_v;
        }

        return bridges;
    }

    vector<vector<int>> get_components(){
        vector<char> is_bridge(m, 0);
        for (auto& b : get_bridges()) is_bridge[b.id] = 1;

        count = 0;
        vector<int> dfs_stack;
        fill(comp.begin(), comp.end(), -1);
        for (int root = 0; root < n; root++){
            if (comp[root] != -1) continue;
            comp[root] = count;
            dfs_stack.push_back(root);
            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                dfs_stack.pop_back();
                for (auto edge : adj[u]){
                    int v = edge.first, id = edge.second;
                    if (comp[v] == -1 && !is_bridge[id]) comp[v] = count, dfs_stack.push_back(v);
                }
            }
            count++;
        }

        vector<vector<int>> groups(count);
        for (int v = 0; v < n; v++) groups[comp[v]].push_back(v);
        return groups;
    }

    vector<vector<int>> get_bridge_tree(){
        get_components();
        vector<vector<int>> tree(count);
        for (auto& b : get_bridges()) tree[comp[b.u]].push_back(comp[b.v]), tree[comp[b.v]].push_back(comp[b.u]);
        return tree;
    }
};

int main(){
    int n, m;
    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        BridgeGraph g(n);
        while (m--){
            int a, b;
            scanf("%d %d", &a, &b);
            g.add_edge(a, b);
        }

        vector<pair<int, int>> weak;
        for (auto& b : g.get_bridges()) weak.push_back({min(b.u, b.v), max(b.u, b.v)});
        sort(weak.begin(), weak.end());

        printf("%d", (int)weak.size());
        for (auto w : weak) printf(" %d %d", w.first, w.second);
        puts("");
    }
    return 0;
}
