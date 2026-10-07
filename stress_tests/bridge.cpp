#include "common.h"

#define main library_main
#include "../code_library/bridge.cpp"
#undef main

/// Component label of every node using all edges except the skipped ones
vector<int> labels(int n, const vector<pair<int, int>>& edges, const set<pair<int, int>>& skipped){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (auto [u, v] : edges){
        if (!skipped.count({min(u, v), max(u, v)})) parent[find(u)] = find(v);
    }
    vector<int> res(n);
    for (int i = 0; i < n; i++) res[i] = find(i);
    return res;
}

/// Simple graphs: no self-loops, no parallel edges
vector<pair<int, int>> random_graph(int n){
    set<pair<int, int>> edges;
    auto add = [&](int u, int v){ if (u != v) edges.insert({min(u, v), max(u, v)}); };
    int density = stress::rand_int(0, n > 12 ? 3 : 100), shape = stress::rand_int(0, 2);

    if (shape == 2 && n >= 2){
        /// Short cycles joined by single edges, a mix of bridges and 2-edge-connected blobs
        for (int c = stress::rand_int(1, n); c; c--){
            vector<int> cycle(stress::rand_int(2, min(n, 5)));
            for (auto& x : cycle) x = stress::rand_int(0, n - 1);
            for (size_t i = 0; i < cycle.size(); i++) add(cycle[i], cycle[(i + 1) % cycle.size()]);
        }
    }
    for (int u = 0; u < n && shape != 2; u++){
        for (int v = u + 1; v < n; v++){
            if (stress::rand_int(0, 99) < density || (shape == 1 && v == u + 1)) add(u, v);
        }
    }

    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    vector<pair<int, int>> res;
    for (auto [u, v] : edges) res.push_back(stress::rand_int(0, 1) ? make_pair(label[u], label[v]) : make_pair(label[v], label[u]));
    shuffle(res.begin(), res.end(), stress::rng());
    return res;
}

int main(){
    /// One graph reused, constructing 100010 adjacency vectors per instance would dominate the run
    static Graph* g = new Graph();

    for (long long it = 0; it < stress::scaled(2500); it++){
        int n = stress::rand_int(1, it % 40 ? 12 : 250);
        auto edges = random_graph(n);
        g->n = n;
        for (int i = 0; i < n; i++) g->adj[i].clear();
        for (auto [u, v] : edges) g->add_edge(u, v);

        auto base = labels(n, edges, {});
        set<pair<int, int>> expected;
        for (auto [u, v] : edges){
            auto cut = labels(n, edges, {{min(u, v), max(u, v)}});
            if (cut[u] != cut[v]) expected.insert({min(u, v), max(u, v)});
        }

        auto bridges = g->get_bridges();
        set<pair<int, int>> found;
        auto split = labels(n, edges, expected);
        for (auto& b : bridges){
            assert(found.insert({min(b.u, b.v), max(b.u, b.v)}).second);
            auto cut = labels(n, edges, {{min(b.u, b.v), max(b.u, b.v)}});
            assert(b.cnt_u == count(cut.begin(), cut.end(), cut[b.u]));
            assert(b.cnt_v == count(cut.begin(), cut.end(), cut[b.v]));
        }
        assert(found == expected);

        /// Bridge tree: nodes share a label exactly when they stay connected without the bridges
        auto tree = g->get_bridge_tree();
        assert(tree.size() == bridges.size());
        map<int, int> to_num, to_split;
        for (int x = 0; x < n; x++){
            assert(to_num.emplace(split[x], g->num[x]).first->second == g->num[x]);
            assert(to_split.emplace(g->num[x], split[x]).first->second == split[x]);
        }
        assert(to_split.empty() || (to_split.begin()->first == 0 && to_split.rbegin()->first == (int)to_split.size() - 1));
        for (size_t i = 0; i < bridges.size(); i++) assert(tree[i] == Pair(g->num[bridges[i].u], g->num[bridges[i].v]));
    }
    return 0;
}
