/***
 *
 * Bridges and the Bridge Tree
 * Bridges, 2-edge-connected components and a strong orientation of an undirected graph, multi-edges and loops allowed
 *
 * Complexity: O(n + m)
 *
 * BridgeGraph g(n): nodes are 0 to n - 1
 * g.add_edge(u, v): adds an edge and returns its id, ids count up from 0 in add_edge order
 * g.get_bridge_tree(): adjacency lists over the g.count components, one tree edge per bridge, a forest
 * g.get_bridges(): each edge whose removal adds a component, with endpoints u, v, id and node counts cnt_u, cnt_v
 * g.get_components(): the 2-edge-connected components, the groups left after removing every bridge, each ascending
 * g.orient(): head[id], the endpoint edge id points to, so every component is strongly connected (details below)
 * g.comp[v]: component of v in [0, g.count) by smallest vertex, set by get_bridge_tree, get_components and orient
 *
 * Notes:
 *   - A doubled edge is never a bridge
 *   - orient points every bridge toward the largest component of its connected part
 *   - That maximizes the fewest vertices any single vertex reaches over all orientations
 *   - A connected graph has a strongly connected orientation iff it has no bridge (Robbins), and orient finds it
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Bridge{
    int u, v;
    int cnt_u, cnt_v; /// node counts on the side of u and of v once the bridge is removed
    int id;
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

    vector<vector<int>> get_bridge_tree(){
        get_components();

        vector<vector<int>> tree(count);
        for (int u = 0; u < n; u++){
            for (auto [v, id] : adj[u]){
                if (comp[u] != comp[v]) tree[comp[u]].push_back(comp[v]);
            }
        }

        return tree;
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
                    auto [v, id] = adj[u][edge_pos[u]++];
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
                if (low[u] > discover[p]) bridges.push_back({p, u, 0, dt - discover[u] + 1, parent_edge[u]});
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
                for (auto [v, id] : adj[u]){
                    if (comp[v] == -1 && !is_bridge[id]) comp[v] = count, dfs_stack.push_back(v);
                }
            }
            count++;
        }

        vector<vector<int>> groups(count);
        for (int v = 0; v < n; v++) groups[comp[v]].push_back(v);
        return groups;
    }

    /// DFS from a vertex of the largest component: tree edges point down, other edges up, bridges back toward the root
    vector<int> orient(){
        vector<vector<int>> groups = get_components();
        vector<char> is_bridge(m, 0), seen(m, 0), visited(n, 0);
        for (auto& b : get_bridges()) is_bridge[b.id] = 1;

        vector<int> head(m), order(n), edge_pos(n, 0), dfs_stack;
        iota(order.begin(), order.end(), 0);
        auto larger = [&](int a, int b){ return groups[comp[a]].size() > groups[comp[b]].size(); };
        stable_sort(order.begin(), order.end(), larger);

        for (int root : order){
            if (visited[root]) continue;
            visited[root] = 1;
            dfs_stack.push_back(root);
            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                if (edge_pos[u] == (int)adj[u].size()){
                    dfs_stack.pop_back();
                    continue;
                }

                auto [v, id] = adj[u][edge_pos[u]++];
                if (seen[id]) continue;
                seen[id] = 1;
                head[id] = !visited[v] && is_bridge[id] ? u : v;
                if (!visited[v]) visited[v] = 1, dfs_stack.push_back(v);
            }
        }

        return head;
    }
};

int main(){
    /// Triangles 0-1-2, 2-3-4 and 5-6-7, the bridge 4-5 between them, and the separate edge 8-9
    BridgeGraph g(10);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);
    g.add_edge(5, 6);
    g.add_edge(6, 7);
    g.add_edge(7, 5);
    g.add_edge(2, 3);
    g.add_edge(3, 4);
    g.add_edge(2, 4);
    g.add_edge(4, 5);
    g.add_edge(8, 9);

    auto bridges = g.get_bridges();
    assert(bridges.size() == 2u);
    assert(bridges[0].id == 9);
    assert(bridges[0].cnt_u == 5);
    assert(bridges[0].cnt_v == 3);
    assert(bridges[1].id == 10);

    assert((g.get_components() == vector<vector<int>>{{0, 1, 2, 3, 4}, {5, 6, 7}, {8}, {9}}));
    assert((g.comp == vector<int>{0, 0, 0, 0, 0, 1, 1, 1, 2, 3}));

    auto tree = g.get_bridge_tree();
    assert(g.count == 4);
    assert((tree == vector<vector<int>>{{1}, {0}, {3}, {2}}));

    /// Closing the outer cycle 0-8-9-7 leaves no bridge, so every node joins component 0
    g.add_edge(8, 0);
    g.add_edge(9, 7);
    assert((g.get_bridge_tree() == vector<vector<int>>{{}}));
    assert(g.comp == vector<int>(10, 0));

    /// Orienting: the two triangles and the 2-3-4 triangle become directed cycles, bridge 4-5 points to {0, 1, 2, 3, 4}
    BridgeGraph roads(8);
    vector<pair<int, int>> streets = {{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 2}, {4, 5}, {5, 6}, {6, 7}, {7, 5}};
    for (auto [u, v] : streets) roads.add_edge(u, v);
    vector<int> head = roads.orient();
    assert(head[6] == 4);
    assert(head[0] == 1);

    /// The doubled edge 0-1 is no bridge and the self loop on 2 changes nothing
    BridgeGraph multi(4);
    multi.add_edge(0, 1);
    multi.add_edge(0, 1);
    multi.add_edge(1, 2);
    multi.add_edge(2, 2);
    multi.add_edge(2, 3);

    auto multi_bridges = multi.get_bridges();
    assert(multi_bridges.size() == 2u);
    assert(multi_bridges[0].id == 4);
    assert(multi_bridges[1].id == 2);

    return 0;
}
