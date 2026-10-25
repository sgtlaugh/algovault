/***
 *
 * Articulation Points, Biconnected Components and the Block-Cut Tree
 * Splits an undirected graph into its blocks (vertex-biconnected components), parallel edges and self loops allowed
 *
 * Complexity: O(n + m)
 *
 * Biconnected g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v): adds an edge and returns its id, ids count up from 0 in add_edge order
 * g.run(): fills blocks, cuts and is_cut, rerun it after more add_edge calls
 * g.block_cut_tree(): adjacency lists over n + blocks nodes, node v < n is vertex v, node n + b is block b, after run()
 * g.blocks: edge ids of each block, two edges share a block iff they lie on a common simple cycle
 * g.cuts: sorted articulation points, the vertices whose removal increases the number of connected components
 * g.is_cut[v]: whether v is an articulation point
 *
 * Notes:
 *   - Every block links to each vertex it touches, so cut vertices have degree >= 2 and the rest are leaves or isolated
 *   - The block-cut tree is a forest with one tree per connected component
 *   - The vertices separating u from w are exactly the cut vertices strictly inside the tree path from u to w
 *   - A bridge is a block of its own, a doubled bridge is a 2-edge block, a self loop belongs to no block
 *   - A vertex lies in several blocks iff it is an articulation point
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Biconnected{
    int n;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge id)
    vector<pair<int, int>> edges;
    vector<vector<int>> blocks;
    vector<int> cuts;
    vector<bool> is_cut;

    Biconnected(int n) : n(n), adj(n) {}

    int add_edge(int u, int v){
        int id = edges.size();
        edges.push_back({u, v});
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        return id;
    }

    void run(){
        vector<int> discover(n, 0), low(n), parent_edge(n, -1), next_edge(n, 0), path, edge_stack;

        is_cut.assign(n, false);
        blocks.clear(), cuts.clear();

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

        for (int u = 0; u < n; u++){
            if (is_cut[u]) cuts.push_back(u);
        }
    }

    vector<vector<int>> block_cut_tree(){
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
};

int main(){
    /// Triangles 0-1-2 and 2-3-4, bridges 4-5 and 5-6, the doubled edge 6-7 and the separate edge 8-9
    Biconnected g(10);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 0);
    g.add_edge(5, 6);
    g.add_edge(6, 7);
    g.add_edge(7, 6);
    g.add_edge(2, 3);
    g.add_edge(3, 4);
    g.add_edge(2, 4);
    g.add_edge(4, 5);
    g.add_edge(8, 9);
    g.run();
    assert((g.cuts == vector<int>{2, 4, 5, 6}));
    assert(g.is_cut[4]);
    assert(!g.is_cut[3]);

    auto sorted_blocks = [](vector<vector<int>> blocks){
        for (auto& block: blocks) sort(block.begin(), block.end());
        sort(blocks.begin(), blocks.end());
        return blocks;
    };
    assert((sorted_blocks(g.blocks) == vector<vector<int>>{{0, 1, 2}, {3}, {4, 5}, {6, 7, 8}, {9}, {10}}));

    /// Block node 10 + b links to the vertices of block b, cut vertices sit in two blocks
    auto tree = g.block_cut_tree();
    assert(tree.size() == 16u);
    map<vector<int>, vector<int>> block_vertices;
    for (int b = 0; b < 6; b++){
        auto block = sorted_blocks({g.blocks[b]})[0];
        auto vertices = tree[10 + b];
        sort(vertices.begin(), vertices.end());
        block_vertices[block] = vertices;
    }
    assert((block_vertices == map<vector<int>, vector<int>>{
        {{0, 1, 2}, {0, 1, 2}}, {{3}, {5, 6}}, {{4, 5}, {6, 7}}, {{6, 7, 8}, {2, 3, 4}}, {{9}, {4, 5}}, {{10}, {8, 9}}
    }));
    vector<int> degrees;
    for (int v = 0; v < 10; v++) degrees.push_back(tree[v].size());
    assert((degrees == vector<int>{1, 1, 2, 1, 2, 2, 2, 1, 1, 1}));

    /// A doubled bridge is a 2-edge block, self loops join no block, so 3 and 4 stay isolated in the tree
    Biconnected multi(5);
    multi.add_edge(0, 1);
    multi.add_edge(1, 0);
    multi.add_edge(1, 1);
    multi.add_edge(1, 2);
    multi.add_edge(3, 3);
    multi.run();
    assert((multi.cuts == vector<int>{1}));
    assert((sorted_blocks(multi.blocks) == vector<vector<int>>{{0, 1}, {3}}));

    auto multi_tree = multi.block_cut_tree();
    assert(multi_tree.size() == 7u);
    assert(multi_tree[1].size() == 2u);
    assert(multi_tree[3].empty());
    assert(multi_tree[4].empty());

    return 0;
}
