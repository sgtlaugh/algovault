/***
 *
 * Cactus Graph
 * Recognizes edge and vertex cacti, lists their cycles and builds the cactus tree
 *
 * Complexity: O(n + m)
 *
 * Edge cactus: every edge lies on at most one simple cycle
 * Vertex cactus: every vertex lies on at most one simple cycle, so every vertex cactus is an edge cactus
 * A disconnected graph qualifies when every component does, check connectivity separately if needed
 * Parallel edges and self loops are allowed: a doubled edge is a cycle of length 2, a self loop one of length 1
 *
 * Cactus g(n); g.add_edge(u, v); bool ok = g.build();
 * build() returns g.edge_cactus and sets g.vertex_cactus, the fields below are valid only when it returns true
 * g.cycles[i]: vertices around cycle i, starting at the one closest to its DFS root
 * g.cycle_edges[i][j]: index of the edge joining cycles[i][j] and cycles[i][(j + 1) % len]
 * g.cycle_of[e]: cycle holding edge e, -1 when e lies on no cycle (a bridge)
 * g.tree: forest on n + cycles.size() nodes, node n + i joined to every vertex of cycle i, each bridge kept as is
 *
 * The DFS is iterative, so long paths and cycles cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Cactus{
    int n;
    bool edge_cactus = false, vertex_cactus = false;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge index)
    vector<pair<int, int>> edges;
    vector<vector<int>> cycles, cycle_edges, tree;
    vector<int> cycle_of;

    Cactus(int n) : n(n), adj(n) {}

    int add_edge(int u, int v){
        adj[u].push_back({v, (int)edges.size()});
        adj[v].push_back({u, (int)edges.size()});
        edges.push_back({u, v});
        return edges.size() - 1;
    }

    bool build(){
        edge_cactus = vertex_cactus = false;
        cycles.clear(), cycle_edges.clear(), tree.clear();
        cycle_of.assign(edges.size(), -1);
        if (!find_cycles()) return false;

        edge_cactus = vertex_cactus = true;
        vector<int> on_cycles(n, 0);
        for (auto& cycle : cycles){
            for (int x : cycle) vertex_cactus &= ++on_cycles[x] <= 1;
        }

        tree.assign(n + cycles.size(), {});
        for (int i = 0; i < (int)cycles.size(); i++){
            for (int x : cycles[i]) tree[x].push_back(n + i), tree[n + i].push_back(x);
        }
        for (int e = 0; e < (int)edges.size(); e++){
            auto [u, v] = edges[e];
            if (cycle_of[e] == -1) tree[u].push_back(v), tree[v].push_back(u);
        }
        return true;
    }

private:
    /// Every non-tree edge closes exactly one cycle with tree edges, a tree edge claimed twice breaks the edge cactus
    bool find_cycles(){
        vector<int> parent_edge(n, -1), it(n, 0), stack;
        vector<char> state(n, 0), used(edges.size(), 0);

        for (int s = 0; s < n; s++){
            if (state[s]) continue;
            state[s] = 1;
            stack.push_back(s);
            while (!stack.empty()){
                int u = stack.back();
                if (it[u] == (int)adj[u].size()){
                    state[u] = 2;
                    stack.pop_back();
                    continue;
                }

                auto [v, id] = adj[u][it[u]++];
                if (used[id]) continue;
                used[id] = 1;
                if (!state[v]){
                    parent_edge[v] = id, state[v] = 1;
                    stack.push_back(v);
                }
                /// A finished v has already used all its edges, so an unused edge always leads up to an ancestor
                else if (!close_cycle(u, v, id, parent_edge)) return false;
            }
        }
        return true;
    }

    bool close_cycle(int u, int v, int id, const vector<int>& parent_edge){
        int c = cycles.size();
        vector<int> cycle = {u}, path;
        for (int x = u; x != v; ){
            int e = parent_edge[x];
            if (cycle_of[e] != -1) return false;
            cycle_of[e] = c;
            path.push_back(e);
            x = edges[e].first ^ edges[e].second ^ x;
            cycle.push_back(x);
        }

        reverse(cycle.begin(), cycle.end());
        reverse(path.begin(), path.end());
        cycle_of[id] = c;
        path.push_back(id);
        cycles.push_back(cycle);
        cycle_edges.push_back(path);
        return true;
    }
};

int main(){
    /***
     * Triangle 0-1-2, bridge 2-3, square 3-4-5-6
    ***/
    Cactus g(7);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 3}}) g.add_edge(u, v);
    assert(g.build());
    assert(g.vertex_cactus);
    assert((g.cycles == vector<vector<int>>{{0, 1, 2}, {3, 4, 5, 6}}));
    assert((g.cycle_of == vector<int>{0, 0, 0, -1, 1, 1, 1, 1}));  /// edge 3 is the bridge 2-3
    sort(g.tree[3].begin(), g.tree[3].end());
    assert((g.tree[3] == vector<int>{2, 8}));                    /// the bridge and the square's node 7 + 1

    /// Two triangles sharing vertex 0
    Cactus bowtie(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {3, 4}, {4, 0}}) bowtie.add_edge(u, v);
    assert(bowtie.build());
    assert(!bowtie.vertex_cactus);

    /// A square with the diagonal 0-2 puts edge 0-1 on two cycles
    Cactus diamond(4);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}) diamond.add_edge(u, v);
    assert(!diamond.build());
    return 0;
}
