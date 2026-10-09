#include "../common.h"

#define main library_main
#include "../../code_library/graphs/bridge.cpp"
#undef main

/// Component label of every node using all edges except the skipped edge indices
vector<int> labels(int n, const vector<pair<int, int>>& edges, const set<int>& skipped){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (int i = 0; i < (int)edges.size(); i++){
        if (!skipped.count(i)) parent[find(edges[i].first)] = find(edges[i].second);
    }
    vector<int> res(n);
    for (int i = 0; i < n; i++) res[i] = find(i);
    return res;
}

/// Simple graphs, or multigraphs with parallel edges and self loops
vector<pair<int, int>> random_graph(int n, bool multi){
    vector<pair<int, int>> edges;
    set<pair<int, int>> seen;
    auto add = [&](int u, int v){
        if (!multi && (u == v || !seen.insert({min(u, v), max(u, v)}).second)) return;
        edges.push_back({u, v});
    };
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
    if (multi){
        for (int extra = stress::rand_int(0, 4); extra && !edges.empty(); extra--){
            auto e = edges[stress::rand_int(0, edges.size() - 1)];
            add(e.first, e.second);
        }
        if (stress::rand_int(0, 1)){
            int x = stress::rand_int(0, n - 1);
            add(x, x);
        }
    }

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

int main(){
    /// One graph reused, constructing 100010 adjacency vectors per instance would dominate the run
    static Graph* g = new Graph();

    for (long long it = 0; it < stress::scaled(2500); it++){
        int n = stress::rand_int(1, it % 40 ? 12 : 250);
        auto edges = random_graph(n, it % 2);
        g->n = n, g->m = 0;
        for (int i = 0; i < n; i++) g->adj[i].clear();
        for (int i = 0; i < (int)edges.size(); i++) assert(g->add_edge(edges[i].first, edges[i].second) == i);

        set<int> expected;
        for (int i = 0; i < (int)edges.size(); i++){
            auto cut = labels(n, edges, {i});
            if (cut[edges[i].first] != cut[edges[i].second]) expected.insert(i);
        }

        auto bridges = g->get_bridges();
        set<int> found;
        auto split = labels(n, edges, expected);
        for (auto& b : bridges){
            assert(found.insert(b.id).second);
            assert(make_pair(min(b.u, b.v), max(b.u, b.v)) == make_pair(min(edges[b.id].first, edges[b.id].second), max(edges[b.id].first, edges[b.id].second)));
            auto cut = labels(n, edges, {b.id});
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

    /// Long paths: every edge (i, i + 1) is a bridge, a bridge set hashed by u ^ v put half of them in one bucket and went quadratic
    /// Cut every 1000 nodes because the recursive DFS overflows an 8 MB stack under ASan beyond about 15000 deep
    int n = 100000, len = 1000;
    g->n = n, g->m = 0;
    for (int i = 0; i < n; i++) g->adj[i].clear();
    for (int i = 0; i + 1 < n; i++){
        if ((i + 1) % len) g->add_edge(i, i + 1);
    }
    auto start = chrono::steady_clock::now();
    auto tree = g->get_bridge_tree();
    assert((int)tree.size() == n - n / len);
    for (int x = 0; x < n; x++) assert(g->num[x] == x);
    assert(chrono::steady_clock::now() - start < chrono::seconds(2));
    return 0;
}
