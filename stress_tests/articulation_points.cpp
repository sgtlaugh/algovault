#include "common.h"

#define main library_main
#include "../code_library/articulation_points.cpp"
#undef main

int components(int n, const vector<pair<int, int>>& edges, int removed){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };

    int res = n - (removed >= 0);
    for (auto [u, v] : edges){
        if (u == removed || v == removed) continue;
        int a = find(u), b = find(v);
        if (a != b) parent[a] = b, res--;
    }
    return res;
}

/// Simple graphs: no self-loops, no parallel edges
vector<pair<int, int>> random_graph(int n){
    set<pair<int, int>> edges;
    int density = stress::rand_int(0, n > 12 ? 3 : 100), shape = stress::rand_int(0, 3);
    auto add = [&](int u, int v){ if (u != v) edges.insert({min(u, v), max(u, v)}); };

    if (shape == 3 && n >= 2){
        /// Short cycles glued at shared nodes, the bowties where back edge handling decides a cut
        for (int c = stress::rand_int(1, n); c; c--){
            vector<int> cycle(stress::rand_int(2, min(n, 5)));
            for (auto& x : cycle) x = stress::rand_int(0, n - 1);
            for (size_t i = 0; i < cycle.size(); i++) add(cycle[i], cycle[(i + 1) % cycle.size()]);
        }
    }
    for (int u = 0; u < n && shape != 3; u++){
        for (int v = u + 1; v < n; v++){
            bool take = stress::rand_int(0, 99) < density;
            if (shape == 1) take = take && v == u + 1;                     /// paths with gaps, every inner node is a cut
            if (shape == 2) take = v == u + 1 || (take && density < 15);  /// a path with a few chords
            if (take) add(u, v);
        }
    }

    vector<pair<int, int>> res;
    for (auto [u, v] : edges) res.push_back(stress::rand_int(0, 1) ? make_pair(u, v) : make_pair(v, u));
    shuffle(res.begin(), res.end(), stress::rng());

    /// Relabel so structure does not follow node order
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    for (auto& [u, v] : res) u = label[u], v = label[v];
    return res;
}

int main(){
    for (long long it = 0; it < stress::scaled(8000); it++){
        int n = stress::rand_int(1, it % 20 ? 12 : 250);
        auto edges = random_graph(n);

        /// One graph reused, constructing 100010 adjacency vectors per instance would dominate the run
        static Graph* g = new Graph();
        g->n = n;
        for (int i = 0; i < n; i++) g->adj[i].clear();

        /// Query halfway through the edges too, a second call must not keep cuts from the first
        int half = stress::rand_int(0, edges.size());
        for (int i = 0; i < (int)edges.size(); i++){
            if (i == half) g->get_cuts();
            g->add_edge(edges[i].first, edges[i].second);
        }

        vector<int> expected;
        int base = components(n, edges, -1);
        for (int v = 0; v < n; v++){
            if (components(n, edges, v) > base) expected.push_back(v);
        }
        assert(g->get_cuts() == expected);
    }
    return 0;
}
