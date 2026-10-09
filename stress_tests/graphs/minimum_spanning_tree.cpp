#include "../common.h"

#define main library_main
#include "../../code_library/graphs/minimum_spanning_tree.cpp"
#undef main

/// Prim's algorithm on an adjacency matrix, run from every unvisited vertex so it covers each component
pair<long long, int> prim(int n, const vector<array<long long, 3>>& edges){
    const long long INF = LLONG_MAX;
    vector<vector<long long>> w(n, vector<long long>(n, INF));
    for (auto [u, v, c] : edges){
        if (u != v) w[u][v] = w[v][u] = min(w[u][v], c);
    }

    vector<char> used(n, 0);
    vector<long long> best(n, INF);
    long long total = 0;
    int components = 0;

    for (int it = 0; it < n; it++){
        int u = -1;
        for (int v = 0; v < n; v++){
            if (!used[v] && (u == -1 || best[v] < best[u])) u = v;
        }
        if (best[u] == INF) components++;
        else total += best[u];
        used[u] = 1;
        for (int v = 0; v < n; v++){
            if (!used[v] && w[u][v] < best[v]) best[v] = w[u][v];
        }
    }
    return {total, components};
}

void check(int n, const vector<array<long long, 3>>& edges, pair<long long, int> expected){
    MSTResult r = minimum_spanning_tree(n, edges);
    auto [total, components] = expected;
    assert(r.weight == total);
    assert((int)r.chosen.size() == n - components);
    assert(r.connected == (components <= 1));

    /// The chosen edges must form a forest: union-find never sees a cycle
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    };

    long long sum = 0;
    for (int i : r.chosen){
        int a = find(edges[i][0]), b = find(edges[i][1]);
        assert(a != b);
        parent[a] = b, sum += edges[i][2];
    }
    assert(sum == r.weight);
}

int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(0, 25), m = n ? stress::rand_int(0, 3 * n) : 0;
        vector<array<long long, 3>> edges;
        for (int i = 0; i < m; i++){
            long long range = it % 3 == 0 ? 3 : 1000000000000LL;
            edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(-range, range)});
        }
        check(n, edges, prim(n, edges));
    }

    /// A large instance with a known answer: tree edge (parent[i] < i, i) weighs i * K - C, so the heaviest tree edge
    /// on any u-v path ends at max(u, v) and an extra edge weighing at least max(u, v) * K - C can never improve the tree
    /// The last 3 vertices stay isolated, and vertices are relabeled so indices carry no order
    const long long K = 1000, C = 100000000;
    int n = 200000, isolated = 3, core = n - isolated;
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<array<long long, 3>> edges;
    long long total = 0;
    for (int i = 1; i < core; i++){
        int p = i % 3 ? i - 1 : stress::rand_int(0, i - 1);
        edges.push_back({label[p], label[i], i * K - C});
        total += i * K - C;
    }
    for (int i = 0; i < n; i++){
        int u = stress::rand_int(0, core - 1), v = stress::rand_int(0, core - 1);
        edges.push_back({label[u], label[v], max(u, v) * K - C + stress::rand_int(0, i % 2 ? 0 : 50 * K)});
    }

    shuffle(edges.begin(), edges.end(), stress::rng());
    check(n, edges, {total, 1 + isolated});
    return 0;
}
