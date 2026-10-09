#include "../common.h"

#define main library_main
#include "../../code_library/trees/dynamic_diameter.cpp"
#undef main

/// Distances from src by an iterative traversal, valid for any non-negative weights
vector<long long> distances(int n, const vector<array<long long, 3>>& edges, int src){
    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges) adj[u].push_back({(int)v, w}), adj[v].push_back({(int)u, w});

    vector<long long> dist(n, -1);
    vector<int> stack = {src};
    dist[src] = 0;
    while (!stack.empty()){
        int u = stack.back();
        stack.pop_back();
        for (auto [v, w] : adj[u]){
            if (dist[v] != -1) continue;
            dist[v] = dist[u] + w;
            stack.push_back(v);
        }
    }
    return dist;
}

/// Double BFS: the farthest node from any start is a diameter endpoint
long long naive_diameter(int n, const vector<array<long long, 3>>& edges){
    auto d0 = distances(n, edges, 0);
    int far = max_element(d0.begin(), d0.end()) - d0.begin();
    auto d1 = distances(n, edges, far);
    return *max_element(d1.begin(), d1.end());
}

/// Random tree of a given shape with random labels
vector<array<long long, 3>> make_tree(int n, int shape, long long max_w){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<array<long long, 3>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        edges.push_back({label[i], label[p], stress::rand_int(0, max_w)});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Applies random updates and compares with double BFS after every check_every-th one
void check(int n, int shape, int updates, long long max_w, int check_every){
    auto edges = make_tree(n, shape, max_w);
    DynamicDiameter tree(n);
    for (int i = 0; i < n - 1; i++) assert(tree.add_edge(edges[i][0], edges[i][1], edges[i][2]) == i);
    tree.build();
    assert(tree.diameter() == naive_diameter(n, edges));
    if (n == 1) return;

    for (int q = 1; q <= updates; q++){
        int id = stress::rand_int(0, n - 2);
        long long w = stress::rand_int(0, 3) ? stress::rand_int(0, max_w) : 0;
        edges[id][2] = w;
        tree.update(id, w);
        if (q % check_every == 0) assert(tree.diameter() == naive_diameter(n, edges));
    }
    assert(tree.diameter() == naive_diameter(n, edges));
}

int main(){
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 400 ? it % 12 + 1 : stress::rand_int(1, 40);
        check(n, it % 5, 40, it % 3 == 0 ? 3 : 1000000000LL, 1);
    }

    for (long long it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(2, 40);
        check(n, it % 5, 40, 2000000000000000000LL / (n - 1), 1);  /// weight sum at the 2e18 limit
    }

    for (long long it = 0; it < stress::scaled(5); it++){
        check(stress::rand_int(20000, 60000), it % 5, 20000, 1000000000LL, 1000);
    }

    check(200000, 1, 20000, 1000000000LL, 5000);  /// a path, deep enough to overflow an 8 MB stack with a recursive tour

    return 0;
}
