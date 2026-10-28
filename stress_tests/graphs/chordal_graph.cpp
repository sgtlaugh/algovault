#include "../common.h"

#define main library_main
#include "../../code_library/graphs/chordal_graph.cpp"
#undef main

vector<pair<int, int>> shuffled(int n, vector<pair<int, int>> edges){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    for (int i = n - 1; i > 0; i--) swap(label[i], label[stress::rand_int(0, i)]);
    for (auto& [u, v] : edges){
        u = label[u], v = label[v];
        if (stress::rand_int(0, 1)) swap(u, v);
    }
    for (int i = (int)edges.size() - 1; i > 0; i--) swap(edges[i], edges[stress::rand_int(0, i)]);

    return edges;
}

/// Each new vertex joins a clique made of an old vertex and part of the clique that old vertex joined
vector<pair<int, int>> random_chordal(int n, int max_clique){
    vector<vector<int>> joined(n);
    vector<pair<int, int>> edges;
    for (int v = 1; v < n; v++){
        if (stress::rand_int(0, 9) == 0) continue;
        int u = stress::rand_int(0, v - 1);
        joined[v] = {u};
        for (int x : joined[u]){
            if ((int)joined[v].size() < max_clique && stress::rand_int(0, 2)) joined[v].push_back(x);
        }
        for (int x : joined[v]) edges.push_back({x, v});
    }

    return shuffled(n, edges);
}

vector<pair<int, int>> random_graph(int n, int percent){
    vector<pair<int, int>> edges;
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (stress::rand_int(0, 99) < percent) edges.push_back({u, v});
        }
    }

    return edges;
}

/// The library checker's test: each vertex's later neighbors all touch its first later neighbor
void check_order(const ChordalGraph& g, const set<pair<int, int>>& edge_set){
    int n = g.n;
    assert((int)g.order.size() == n);
    vector<int> pos(n, -1);
    for (int i = 0; i < n; i++){
        assert(0 <= g.order[i] && g.order[i] < n && pos[g.order[i]] == -1);
        pos[g.order[i]] = i;
    }

    for (int v = 0; v < n; v++){
        int first = -1;
        for (int w : g.adj[v]){
            if (pos[w] > pos[v] && (first == -1 || pos[w] < pos[first])) first = w;
        }
        for (int w : g.adj[v]){
            if (pos[w] > pos[v] && w != first) assert(edge_set.count(minmax(first, w)));
        }
    }
}

/// Induced: the cycle's consecutive pairs are edges and no other edge has both ends on it
void check_cycle(const ChordalGraph& g, const vector<pair<int, int>>& edges, const set<pair<int, int>>& edge_set){
    int k = g.cycle.size();
    assert(k >= 4);
    vector<char> on_cycle(g.n, 0);
    for (int i = 0; i < k; i++){
        assert(!on_cycle[g.cycle[i]]);
        on_cycle[g.cycle[i]] = 1;
        assert(edge_set.count(minmax(g.cycle[i], g.cycle[(i + 1) % k])));
    }

    int inside = 0;
    for (auto [u, v] : edges) inside += on_cycle[u] && on_cycle[v];
    assert(inside == k);
}

/// Validates whichever certificate run() returns, on any graph
ChordalGraph check(int n, const vector<pair<int, int>>& edges){
    ChordalGraph g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);
    set<pair<int, int>> edge_set;
    for (auto [u, v] : edges) edge_set.insert(minmax(u, v));

    if (g.run()){
        assert(g.cycle.empty());
        check_order(g, edge_set);
    }
    else check_cycle(g, edges, edge_set);

    return g;
}

/// Brute force on n <= 8: chordality from every vertex subset, clique and independence numbers likewise
void check_small(int n, const vector<pair<int, int>>& edges){
    vector<int> adj(n, 0);
    for (auto [u, v] : edges) adj[u] |= 1 << v, adj[v] |= 1 << u;

    bool has_hole = false;
    int clique = 0, independent = 0;
    for (int s = 0; s < (1 << n); s++){
        bool is_clique = true, is_independent = true, all_degree_two = true;
        for (int v = 0; v < n; v++){
            if (!(s >> v & 1)) continue;
            int inside = __builtin_popcount(adj[v] & s);
            is_clique &= inside == __builtin_popcount(s) - 1;
            is_independent &= inside == 0;
            all_degree_two &= inside == 2;
        }
        if (is_clique) clique = max(clique, __builtin_popcount(s));
        if (is_independent) independent = max(independent, __builtin_popcount(s));
        if (!all_degree_two || __builtin_popcount(s) < 4) continue;

        int reach = s & -s;
        for (int step = 0; step < n; step++){
            for (int v = 0; v < n; v++){
                if (reach >> v & 1) reach |= adj[v] & s;
            }
        }
        has_hole |= reach == s;
    }

    ChordalGraph g = check(n, edges);
    assert(g.cycle.empty() == !has_hole);
    if (has_hole) return;

    vector<int> c = g.max_clique(), is = g.max_independent_set();
    assert((int)c.size() == clique && (int)is.size() == independent);
    for (int x : c){
        for (int y : c) assert(x == y || (adj[x] >> y & 1));
    }
    for (int x : is){
        for (int y : is) assert(!(adj[x] >> y & 1));
    }
}

int main(){
    for (int n = 0; n <= 6; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }
        for (int mask = 0; mask < (1 << pairs.size()); mask++){
            vector<pair<int, int>> edges;
            for (int i = 0; i < (int)pairs.size(); i++){
                if (mask >> i & 1) edges.push_back(pairs[i]);
            }
            check_small(n, shuffled(n, edges));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(7, 8);
        if (it % 2) check_small(n, random_graph(n, stress::rand_int(10, 90)));
        else check_small(n, random_chordal(n, stress::rand_int(1, n)));
    }

    /// Medium graphs: chordal ones must be recognized, a few extra edges usually break them
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(10, 300);
        auto edges = random_chordal(n, stress::rand_int(1, 6));
        assert(check(n, edges).cycle.empty());

        set<pair<int, int>> edge_set;
        for (auto [u, v] : edges) edge_set.insert(minmax(u, v));
        for (int extra = stress::rand_int(1, 3); extra > 0; extra--){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            if (u != v && edge_set.insert(minmax(u, v)).second) edges.push_back({u, v});
        }
        check(n, edges);
        check(n, random_graph(min(n, 60), stress::rand_int(1, 30)));
    }

    /// A cycle through all n vertices is its own only certificate, and without one edge it is a chordal path
    int n = 100000;
    vector<pair<int, int>> ring;
    for (int i = 0; i < n; i++) ring.push_back({i, (i + 1) % n});
    assert((int)check(n, shuffled(n, ring)).cycle.size() == n);

    ring.pop_back();
    assert(check(n, shuffled(n, ring)).cycle.empty());
    assert(check(n, random_chordal(n, 20)).cycle.empty());
    return 0;
}
