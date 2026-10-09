#include "../common.h"

#define main library_main
#include "../../code_library/graphs/kruskal_reconstruction_tree.cpp"
#undef main

vector<vector<int>> light_adjacency(int n, const vector<array<long long, 3>>& edges, long long w){
    vector<vector<int>> adj(n);
    for (auto [a, b, c] : edges){
        if (c <= w) adj[a].push_back(b), adj[b].push_back(a);
    }
    return adj;
}

/// Sorted so it compares directly against the krt range
vector<int> reachable(const vector<vector<int>>& adj, int v){
    vector<int> queue = {v};
    vector<bool> seen(adj.size(), false);
    seen[v] = true;
    for (int i = 0; i < (int)queue.size(); i++){
        for (int u : adj[queue[i]]){
            if (!seen[u]) seen[u] = true, queue.push_back(u);
        }
    }

    sort(queue.begin(), queue.end());
    return queue;
}

vector<int> component(const KruskalReconstructionTree<long long>& krt, int v, long long w){
    auto [l, r] = krt.range(v, w);
    assert(0 <= l && l < r && r <= krt.n);
    assert(krt.size(v, w) == r - l);
    vector<int> res(krt.order.begin() + l, krt.order.begin() + r);
    sort(res.begin(), res.end());
    return res;
}

KruskalReconstructionTree<long long> make_krt(int n, const vector<array<long long, 3>>& edges){
    KruskalReconstructionTree<long long> krt(n);
    for (auto [u, v, w] : edges) krt.add_edge(u, v, w);
    krt.build();
    for (int v = 0; v < n; v++) assert(krt.order[krt.pos[v]] == v);
    return krt;
}

void check(int n, const vector<array<long long, 3>>& edges, const vector<pair<int, long long>>& queries){
    auto krt = make_krt(n, edges);
    for (auto [v, w] : queries) assert(component(krt, v, w) == reachable(light_adjacency(n, edges, w), v));
}

/// Every vertex against every threshold at and around each weight, one BFS per component per threshold
void check_all_thresholds(int n, const vector<array<long long, 3>>& edges){
    auto krt = make_krt(n, edges);
    set<long long> thresholds = {LLONG_MIN, LLONG_MAX};
    for (auto& e : edges) thresholds.insert(e[2] - 1), thresholds.insert(e[2]), thresholds.insert(e[2] + 1);

    for (long long w : thresholds){
        auto adj = light_adjacency(n, edges, w);
        vector<vector<int>> expected(n);
        for (int v = 0; v < n; v++){
            if (!expected[v].empty()) continue;
            auto comp = reachable(adj, v);
            for (int u : comp) expected[u] = comp;
        }
        for (int v = 0; v < n; v++) assert(component(krt, v, w) == expected[v]);
    }
}

vector<array<long long, 3>> random_edges(int n, int m, long long lo, long long hi){
    vector<array<long long, 3>> edges;
    for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(lo, hi)});
    return edges;
}

int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it % 12 + 1 : stress::rand_int(1, 30);
        int m = stress::rand_int(0, 2 * n);
        long long spread = it % 3 == 0 ? 2 : it % 3 == 1 ? 10 : 1000000000000000000LL;
        check_all_thresholds(n, random_edges(n, m, -spread, spread));
    }

    for (long long it = 0; it < stress::scaled(20); it++){
        int n = stress::rand_int(2000, 5000);
        auto edges = random_edges(n, stress::rand_int(n / 2, 2 * n), -1000, 1000);
        vector<pair<int, long long>> queries;
        for (int q = 0; q < 50; q++) queries.push_back({stress::rand_int(0, n - 1), stress::rand_int(-1001, 1001)});
        check(n, edges, queries);
    }

    /// a path with increasing weights builds a tree of depth n - 1, the deepest binary lifting has to cover
    int n = 200000;
    vector<array<long long, 3>> path;
    for (int i = 0; i + 1 < n; i++) path.push_back({i, i + 1, i});
    vector<pair<int, long long>> queries = {{0, -1}, {0, 0}, {n - 1, n - 3}, {n - 1, n - 2}, {n / 2, n / 3}, {0, n}};
    for (int q = 0; q < 20; q++) queries.push_back({stress::rand_int(0, n - 1), stress::rand_int(-1, n)});
    check(n, path, queries);

    check(n, random_edges(n, 2 * n, 1, 1000000000), {{0, 1}, {1, 500000000}, {2, 1000000000}, {3, 20000}});

    return 0;
}
