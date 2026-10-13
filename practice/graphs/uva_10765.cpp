/***
 *
 * Articulation Points, Biconnected Components and the Block-Cut Tree
 * Splits an undirected graph into its blocks (vertex-biconnected components), parallel edges and self loops allowed
 *
 * Complexity: O(n + m)
 *
 * Biconnected g(n); int id = g.add_edge(u, v); nodes are 0 to n - 1, edge ids count up from 0 in add_edge order
 * g.get_cuts(): sorted articulation points, the vertices whose removal increases the number of connected components
 * g.get_blocks(): edge ids of each block, two edges share a block iff they lie on a common simple cycle
 *     A bridge is a block of its own, a doubled bridge is a 2-edge block, a self loop belongs to no block
 *     A vertex can lie in several blocks: exactly the articulation points do
 * g.get_block_cut_tree(): adjacency lists over n + blocks nodes, node v < n is vertex v, node n + b is block b
 *     Every block links to each vertex it touches, so cut vertices have degree >= 2 and the rest are leaves or isolated
 *     The result is a forest with one tree per connected component, the vertices separating u from w are
 *     exactly the cut vertices strictly inside the tree path from u to w
 * Each getter reruns g.decompose(), which fills g.is_cut and g.blocks, so edges added between calls are seen
 *
 * The DFS is iterative, so long paths cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Biconnected{
    int n;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge id)
    vector<pair<int, int>> edges;
    vector<bool> is_cut;
    vector<vector<int>> blocks;

    Biconnected(int n) : n(n), adj(n) {}

    int add_edge(int u, int v){
        int id = edges.size();
        edges.push_back({u, v});
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        return id;
    }

    void decompose(){
        vector<int> discover(n, 0), low(n), parent_edge(n, -1), next_edge(n, 0), path, edge_stack;
        is_cut.assign(n, false);
        blocks.clear();

        int disc_t = 0;
        for (int root = 0; root < n; root++){
            if (discover[root]) continue;
            int children = 0;
            discover[root] = low[root] = ++disc_t;
            path.push_back(root);

            while (!path.empty()){
                int u = path.back();
                if (next_edge[u] < (int)adj[u].size()){
                    auto [v, id] = adj[u][next_edge[u]++];
                    /// Skip only the edge we arrived by, a parallel edge back to the parent is a real back edge
                    if (v == u || id == parent_edge[u]) continue;
                    if (!discover[v]){
                        discover[v] = low[v] = ++disc_t;
                        parent_edge[v] = id;
                        edge_stack.push_back(id);
                        path.push_back(v);
                    }
                    else if (discover[v] < discover[u]){
                        low[u] = min(low[u], discover[v]);
                        edge_stack.push_back(id);
                    }
                    continue;
                }

                path.pop_back();
                if (path.empty()) break;
                int p = path.back();
                low[p] = min(low[p], low[u]);
                if (low[u] < discover[p]) continue;

                if (p != root) is_cut[p] = true;
                else if (++children == 2) is_cut[p] = true;

                blocks.emplace_back();
                int id;
                do{
                    id = edge_stack.back();
                    edge_stack.pop_back();
                    blocks.back().push_back(id);
                } while (id != parent_edge[u]);
            }
        }
    }

    vector<vector<int>> get_block_cut_tree(){
        decompose();

        vector<vector<int>> tree(n + blocks.size());
        vector<int> last(n, -1);
        for (int b = 0; b < (int)blocks.size(); b++){
            for (int id: blocks[b]){
                for (int u: {edges[id].first, edges[id].second}){
                    if (last[u] == b) continue;
                    last[u] = b;
                    tree[u].push_back(n + b);
                    tree[n + b].push_back(u);
                }
            }
        }
        return tree;
    }

    vector<vector<int>> get_blocks(){
        decompose();
        return blocks;
    }

    vector<int> get_cuts(){
        decompose();

        vector<int> cuts;
        for (int u = 0; u < n; u++){
            if (is_cut[u]) cuts.push_back(u);
        }
        return cuts;
    }
};

int main(){
    int n, m, x, y;
    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        Biconnected g(n);
        while (scanf("%d %d", &x, &y) == 2 && x != -1) g.add_edge(x, y);

        /// The network is connected, so bombing a station leaves one piece per block through it: its block-cut tree degree
        auto tree = g.get_block_cut_tree();
        vector<pair<int, int>> order(n);
        for (int v = 0; v < n; v++) order[v] = {-(int)tree[v].size(), v};
        sort(order.begin(), order.end());

        for (int i = 0; i < m; i++) printf("%d %d\n", order[i].second, -order[i].first);
        puts("");
    }
    return 0;
}
