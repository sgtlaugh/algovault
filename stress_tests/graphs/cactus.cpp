#include "../common.h"

#define main library_main
#include "../../code_library/graphs/cactus.cpp"
#undef main

/// Every simple cycle as a sorted list of edge indices, found from its smallest vertex
set<vector<int>> all_cycles(int n, const vector<pair<int, int>>& edges){
    vector<vector<pair<int, int>>> adj(n);
    for (int i = 0; i < (int)edges.size(); i++){
        adj[edges[i].first].push_back({edges[i].second, i});
        adj[edges[i].second].push_back({edges[i].first, i});
    }

    set<vector<int>> found;
    vector<char> on_path(n, 0);
    vector<int> path;
    function<void(int, int)> extend = [&](int s, int x){
        for (auto [y, id] : adj[x]){
            if (y == s && find(path.begin(), path.end(), id) == path.end()){
                auto cycle = path;
                cycle.push_back(id);
                sort(cycle.begin(), cycle.end());
                found.insert(cycle);
            }
            else if (y > s && !on_path[y]){
                on_path[y] = 1, path.push_back(id);
                extend(s, y);
                on_path[y] = 0, path.pop_back();
            }
        }
    };
    for (int s = 0; s < n; s++){
        on_path[s] = 1;
        extend(s, s);
        on_path[s] = 0;
    }
    return found;
}

/// (edge cactus, vertex cactus) straight from the definitions: no edge, resp. vertex, on two simple cycles
pair<bool, bool> cactus_flags(int n, const vector<pair<int, int>>& edges, const set<vector<int>>& cycles){
    vector<int> edge_hits(edges.size(), 0), vertex_hits(n, 0);
    for (auto& cycle : cycles){
        set<int> vertices;
        for (int e : cycle) edge_hits[e]++, vertices.insert(edges[e].first), vertices.insert(edges[e].second);
        for (int x : vertices) vertex_hits[x]++;
    }

    auto twice = [](int h){ return h > 1; };
    return {none_of(edge_hits.begin(), edge_hits.end(), twice), none_of(vertex_hits.begin(), vertex_hits.end(), twice)};
}

/// Cactus grown by bridges, cycles, doubled edges and self loops, then a few random extra edges that usually break it
vector<pair<int, int>> random_graph(int n, int extra, bool multi){
    vector<pair<int, int>> edges;
    int built = 1;
    while (built < n){
        int at = stress::rand_int(0, built - 1), kind = stress::rand_int(0, multi ? 3 : 1);
        if (kind == 0 || built + 1 == n) edges.push_back({at, built++});
        else if (kind == 1){
            int len = stress::rand_int(2, min(5, n - built + 1)), prev = at;
            for (int i = 1; i < len; i++) edges.push_back({prev, built}), prev = built++;
            edges.push_back({prev, at});
        }
        else if (kind == 2){
            edges.push_back({at, built}), edges.push_back({built++, at});
        }
        else edges.push_back({at, at});
    }
    for (int i = 0; i < extra; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (multi || u != v) edges.push_back({u, v});
    }

    /// Several components: drop some edges, cutting a cactus never makes it a non-cactus
    for (int drop = stress::rand_int(0, 2); drop && !edges.empty(); drop--) edges.erase(edges.begin() + stress::rand_int(0, edges.size() - 1));

    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    for (auto& [u, v] : edges){
        u = label[u], v = label[v];
        if (stress::rand_int(0, 1)) swap(u, v);
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

int find_root(vector<int>& parent, int x){
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

/// Cycles really are closed walks over their listed edges, cycle_of agrees, and the tree is a forest with the graph's components
void check_structure(const Cactus& g){
    int n = g.n, k = g.cycles.size(), m = g.edges.size(), bridges = 0;
    vector<int> owner(m, -1);
    for (int i = 0; i < k; i++){
        auto& cycle = g.cycles[i];
        auto& ids = g.cycle_edges[i];
        assert(!cycle.empty() && cycle.size() == ids.size());
        assert(set<int>(cycle.begin(), cycle.end()).size() == cycle.size());
        for (size_t j = 0; j < cycle.size(); j++){
            int a = cycle[j], b = cycle[(j + 1) % cycle.size()];
            auto [u, v] = g.edges[ids[j]];
            assert((u == a && v == b) || (u == b && v == a));
            assert(owner[ids[j]] == -1);
            owner[ids[j]] = i;
        }
    }
    assert(owner == g.cycle_of);

    for (int e = 0; e < m; e++) bridges += owner[e] == -1;
    assert((int)g.tree.size() == n + k);
    vector<int> parent(n + k), graph_parent(n);
    iota(parent.begin(), parent.end(), 0);
    iota(graph_parent.begin(), graph_parent.end(), 0);
    long long degree_sum = 0, merges = 0;
    for (int x = 0; x < n + k; x++){
        degree_sum += g.tree[x].size();
        for (int y : g.tree[x]){
            assert(x < n || y < n);
            int a = find_root(parent, x), b = find_root(parent, y);
            if (a != b) parent[a] = b, merges++;
        }
    }

    long long cycle_vertices = 0;
    for (auto& cycle : g.cycles) cycle_vertices += cycle.size();
    assert(degree_sum == 2 * (bridges + cycle_vertices));
    assert(merges * 2 == degree_sum);
    for (auto [u, v] : g.edges) graph_parent[find_root(graph_parent, u)] = find_root(graph_parent, v);
    map<int, int> tree_to_graph, graph_to_tree;
    for (int x = 0; x < n; x++){
        int a = find_root(parent, x), b = find_root(graph_parent, x);
        assert(tree_to_graph.emplace(a, b).first->second == b);
        assert(graph_to_tree.emplace(b, a).first->second == a);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 9);
        auto edges = random_graph(n, stress::rand_int(0, 2), it % 2);
        Cactus g(n);
        for (int i = 0; i < (int)edges.size(); i++) assert(g.add_edge(edges[i].first, edges[i].second) == i);

        auto expected = all_cycles(n, edges);
        auto [edge_ok, vertex_ok] = cactus_flags(n, edges, expected);
        bool built = g.build();
        assert(built == edge_ok && g.edge_cactus == edge_ok && g.vertex_cactus == vertex_ok);
        if (!built) continue;

        set<vector<int>> listed;
        for (auto ids : g.cycle_edges){
            sort(ids.begin(), ids.end());
            listed.insert(ids);
        }
        assert(listed == expected && listed.size() == g.cycles.size());
        check_structure(g);

        /// build() must reset everything, a second call sees the same graph
        auto cycles = g.cycles;
        auto tree = g.tree;
        assert(g.build() && g.cycles == cycles && g.tree == tree);

        /// Reusing the instance after one more edge, which often breaks the cactus, must not leak the previous flags
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        edges.push_back({u, v});
        g.add_edge(u, v);
        auto [edge_after, vertex_after] = cactus_flags(n, edges, all_cycles(n, edges));
        assert(g.build() == edge_after && g.edge_cactus == edge_after && g.vertex_cactus == vertex_after);
    }

    /// One 150000-cycle followed by a chain of triangles: depth ~250000, an 8 MB stack would not survive a recursive DFS
    int ring = 150000, triangles = 50000, n = ring + 2 * triangles;
    Cactus big(n);
    for (int i = 0; i < ring; i++) big.add_edge(i, (i + 1) % ring);
    for (int t = 0, at = ring - 1; t < triangles; t++){
        int a = ring + 2 * t, b = a + 1;
        big.add_edge(at, a), big.add_edge(a, b), big.add_edge(b, at);
        at = b;
    }

    auto start = chrono::steady_clock::now();
    assert(big.build() && !big.vertex_cactus);
    assert((int)big.cycles.size() == 1 + triangles);
    assert(chrono::steady_clock::now() - start < chrono::seconds(2));
    check_structure(big);

    big.add_edge(0, ring / 2);
    assert(!big.build() && !big.edge_cactus && !big.vertex_cactus);
    return 0;
}
