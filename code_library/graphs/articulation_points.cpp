/***
 *
 * Articulation Points, Biconnected Components and the Block-Cut Tree
 * Splits an undirected graph into its blocks (vertex-biconnected components), parallel edges and self loops allowed
 *
 * Complexity: O(n + m)
 *
 * Graph g(n); int id = g.add_edge(u, v); nodes are 0 to n - 1, edge ids count up from 0 in add_edge order
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

struct Graph{
    int n;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge id)
    vector<pair<int, int>> edges;
    vector<bool> is_cut;
    vector<vector<int>> blocks;

    Graph(int n) : n(n), adj(n) {}

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

    auto sorted_blocks = [](vector<vector<int>> blocks){
        for (auto& block: blocks) sort(block.begin(), block.end());
        sort(blocks.begin(), blocks.end());
        return blocks;
    };

    auto blocks = graph.get_blocks();
    assert((sorted_blocks(blocks) == vector<vector<int>>{{0, 1, 2}, {3}, {4, 5}, {6, 7, 8}, {9}, {10}}));

    auto tree = graph.get_block_cut_tree();
    assert(tree.size() == 16u);
    map<vector<int>, vector<int>> block_vertices;
    for (int b = 0; b < 6; b++){
        auto block = sorted_blocks({blocks[b]})[0];
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

    auto star = Graph(4);
    star.add_edge(0, 1);
    star.add_edge(0, 2);
    star.add_edge(0, 3);
    assert((star.get_cuts() == vector<int>{0}));
    assert((sorted_blocks(star.get_blocks()) == vector<vector<int>>{{0}, {1}, {2}}));

    auto cycle = Graph(4);
    for (int v = 0; v < 4; v++) cycle.add_edge(v, (v + 1) % 4);
    assert(cycle.get_cuts().empty());
    assert((sorted_blocks(cycle.get_blocks()) == vector<vector<int>>{{0, 1, 2, 3}}));

    auto multi = Graph(5);
    multi.add_edge(0, 1);
    multi.add_edge(1, 0);
    multi.add_edge(1, 1);
    multi.add_edge(1, 2);
    multi.add_edge(3, 3);
    assert((multi.get_cuts() == vector<int>{1}));
    assert((sorted_blocks(multi.get_blocks()) == vector<vector<int>>{{0, 1}, {3}}));
    auto multi_tree = multi.get_block_cut_tree();
    assert(multi_tree.size() == 7u);
    assert(multi_tree[3].empty() && multi_tree[4].empty() && multi_tree[1].size() == 2u);

    auto empty = Graph(0);
    assert(empty.get_cuts().empty() && empty.get_blocks().empty() && empty.get_block_cut_tree().empty());

    return 0;
}
