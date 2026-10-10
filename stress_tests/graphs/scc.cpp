#include "../common.h"

#define main library_main
#include "../../code_library/graphs/scc.cpp"
#undef main

/// Same component iff mutually reachable, components numbered along edges, condensation is the exact DAG
void check(int n, const vector<pair<int, int>>& edges){
    /// Run once halfway through the edges too, so the final run must not keep components from the first
    SCC g(n);
    for (int i = 0; i < (int)edges.size(); i++){
        if (i == (int)edges.size() / 2) g.run();
        g.add_edge(edges[i].first, edges[i].second);
    }
    g.run();

    vector<vector<int>> adj(n);
    for (auto [u, v] : edges) adj[u].push_back(v);
    vector<vector<char>> reach(n, vector<char>(n, 0));
    for (int s = 0; s < n; s++){
        vector<int> queue = {s};
        reach[s][s] = 1;
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : adj[queue[i]]){
                if (!reach[s][v]) reach[s][v] = 1, queue.push_back(v);
            }
        }
    }

    set<int> ids;
    for (int u = 0; u < n; u++){
        assert(0 <= g.comp[u] && g.comp[u] < g.count);
        ids.insert(g.comp[u]);
        for (int v = 0; v < n; v++) assert((g.comp[u] == g.comp[v]) == (reach[u][v] && reach[v][u]));
    }
    assert((int)ids.size() == g.count);
    for (auto [u, v] : edges) assert(g.comp[u] <= g.comp[v]);

    set<pair<int, int>> expected;
    for (auto [u, v] : edges){
        if (g.comp[u] != g.comp[v]) expected.insert({g.comp[u], g.comp[v]});
    }

    set<pair<int, int>> got;
    auto dag = g.condensation();
    for (int c = 0; c < g.count; c++){
        assert(is_sorted(dag[c].begin(), dag[c].end()));
        for (int d : dag[c]) assert(got.insert({c, d}).second);
    }
    assert(got == expected);
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(0, 25), m = n ? stress::rand_int(0, 3 * n) : 0;
        vector<pair<int, int>> edges;
        for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        check(n, edges);
    }

    /// One long cycle and one long path, deep enough to overflow a recursive DFS on an 8 MB stack
    int n = 300000;
    SCC cycle(n), path(n);
    for (int i = 0; i < n; i++) cycle.add_edge(i, (i + 1) % n);
    for (int i = 0; i + 1 < n; i++) path.add_edge(i, i + 1);

    cycle.run(), path.run();
    assert(cycle.count == 1 && path.count == n);
    for (int i = 0; i < n; i++) assert(path.comp[i] == i);
    return 0;
}
