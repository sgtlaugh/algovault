#include "../common.h"

#define main library_main
#include "../../code_library/graphs/reachability.cpp"
#undef main

vector<vector<char>> brute_reach(int n, const vector<pair<int, int>>& edges){
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges) adj[u].push_back(v);

    vector<vector<char>> seen(n, vector<char>(n, 0));
    for (int s = 0; s < n; s++){
        vector<int> queue = {s};
        seen[s][s] = 1;
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : adj[queue[i]]){
                if (!seen[s][v]) seen[s][v] = 1, queue.push_back(v);
            }
        }
    }

    return seen;
}

/// shape 0: any edges, 1: a DAG under a random labeling, 2: a few big cycles chained together
vector<pair<int, int>> random_edges(int n, int m, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int i = 0; i < m && n > 0; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (shape >= 1 && u > v) swap(u, v);
        edges.push_back({label[u], label[v]});
    }
    if (shape == 2){
        for (int i = 0; i + 1 < n; i++){
            if (stress::rand_int(0, 7)) edges.push_back({label[i + 1], label[i]});
        }
    }

    return edges;
}

void check(int n, const vector<pair<int, int>>& edges, bool with_closure){
    Reachability g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);
    vector<vector<char>> seen = brute_reach(n, edges);

    vector<int> counts = g.counts();
    for (int u = 0; u < n; u++) assert(counts[u] == count(seen[u].begin(), seen[u].end(), 1));

    vector<pair<int, int>> queries;
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++) queries.push_back({u, v});
    }
    shuffle(queries.begin(), queries.end(), stress::rng());
    vector<char> got = g.reaches(queries);
    for (int i = 0; i < (int)queries.size(); i++) assert(got[i] == seen[queries[i].first][queries[i].second]);

    if (!with_closure) return;
    vector<vector<uint64_t>> reach = g.closure();
    assert((int)reach.size() == n);
    for (int u = 0; u < n; u++){
        assert((int)reach[u].size() == (n + 63) / 64);
        for (int v = 0; v < (int)reach[u].size() * 64; v++){
            bool bit = reach[u][v >> 6] >> (v & 63) & 1;
            assert(bit == (v < n && seen[u][v]));
        }
    }
}

int main(){
    Reachability empty(0);
    assert(empty.counts().empty() && empty.closure().empty() && empty.reaches({}).empty());

    /// Every graph on up to 3 vertices, self loops included
    for (int n = 1; n <= 3; n++){
        for (int mask = 0; mask < 1 << (n * n); mask++){
            vector<pair<int, int>> edges;
            for (int i = 0; i < n * n; i++){
                if (mask >> i & 1) edges.push_back({i / n, i % n});
            }
            check(n, edges, true);
        }
    }

    for (int it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, 140), m = stress::rand_int(0, 3 * n);
        check(n, random_edges(n, m, it % 3), true);
    }

    /// Big graphs cross many 64-vertex chunks, the closure checked on a sample of rows
    for (int it = 0; it < stress::scaled(4); it++){
        int n = stress::rand_int(1500, 2500), m = stress::rand_int(n, 4 * n);
        vector<pair<int, int>> edges = random_edges(n, m, it % 3);
        Reachability g(n);
        for (auto [u, v] : edges) g.add_edge(u, v);
        vector<vector<char>> seen = brute_reach(n, edges);

        vector<int> counts = g.counts();
        vector<vector<uint64_t>> reach = g.closure();
        vector<pair<int, int>> queries;
        for (int i = 0; i < 200000; i++) queries.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        vector<char> got = g.reaches(queries);

        for (int u = 0; u < n; u++) assert(counts[u] == count(seen[u].begin(), seen[u].end(), 1));
        for (int i = 0; i < (int)queries.size(); i++){
            auto [u, v] = queries[i];
            assert(got[i] == seen[u][v]);
            assert((reach[u][v >> 6] >> (v & 63) & 1) == (uint64_t)seen[u][v]);
        }
    }

    return 0;
}
