#include "../common.h"

#define main library_main
#include "../../code_library/graphs/absolute_center.cpp"
#undef main

const long long NONE = LLONG_MAX / 3;

vector<vector<long long>> floyd(int n, const vector<array<long long, 3>>& edges){
    vector<vector<long long>> d(n, vector<long long>(n, NONE));
    for (int v = 0; v < n; v++) d[v][v] = 0;
    for (auto [u, v, w] : edges){
        d[u][v] = min(d[u][v], w);
        d[v][u] = min(d[v][u], w);
    }
    for (int k = 0; k < n; k++){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
        }
    }

    return d;
}

/// Largest shortest distance between two vertices, NONE when disconnected
long long diameter(int n, const vector<array<long long, 3>>& edges){
    vector<vector<long long>> d = floyd(n, edges);
    long long res = 0;
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++) res = max(res, d[u][v]);
    }

    return res;
}

/// Twice the farthest vertex distance from the point at distance offset / 2 from u on edge (u, v, w)
long long farthest(const vector<vector<long long>>& d, array<long long, 3> e, long long offset){
    auto [u, v, w] = e;
    long long res = 0;
    for (int x = 0; x < (int)d.size(); x++) res = max(res, min(2 * d[u][x] + offset, 2 * d[v][x] + 2 * w - offset));
    return res;
}

/// Every offset where a rising tent side meets a falling one, plus both ends: the envelope's minima are among them
long long brute_center(const vector<vector<long long>>& d, const vector<array<long long, 3>>& edges, bool every_offset){
    int n = d.size();
    long long best = n == 1 ? 0 : NONE;
    for (auto e : edges){
        auto [u, v, w] = e;
        vector<long long> offsets = {0, 2 * w};
        if (every_offset){
            for (long long t = 1; t < 2 * w; t++) offsets.push_back(t);
        }
        else{
            for (int x = 0; x < n; x++){
                for (int y = 0; y < n; y++){
                    long long t = d[v][y] + w - d[u][x];
                    if (0 <= t && t <= 2 * w) offsets.push_back(t);
                }
            }
        }
        for (long long t : offsets) best = min(best, farthest(d, e, t));
    }

    return best;
}

/// Asserts ids are n - 1 distinct edges spanning all vertices, and returns the tree's diameter
long long tree_diameter(int n, const vector<array<long long, 3>>& edges, const vector<int>& ids){
    assert((int)ids.size() == n - 1);
    vector<bool> seen(edges.size(), false);
    vector<array<long long, 3>> chosen;
    for (int id : ids){
        assert(0 <= id && id < (int)edges.size() && !seen[id]);
        seen[id] = true;
        chosen.push_back(edges[id]);
    }

    long long res = diameter(n, chosen);
    assert(res < NONE);
    return res;
}

/// Smallest diameter over every n - 1 edge subset that spans: needs no center theory at all
long long brute_mdst(int n, const vector<array<long long, 3>>& edges){
    int m = edges.size();
    long long best = NONE;
    for (int mask = 0; mask < (1 << m); mask++){
        if (__builtin_popcount(mask) != n - 1) continue;
        vector<array<long long, 3>> chosen;
        for (int i = 0; i < m; i++){
            if (mask >> i & 1) chosen.push_back(edges[i]);
        }
        best = min(best, diameter(n, chosen));
    }

    return best;
}

/// A random spanning tree first keeps the graph connected, then extra edges may repeat or loop
vector<array<long long, 3>> random_graph(int n, int extra, long long max_w){
    vector<array<long long, 3>> edges;
    for (int v = 1; v < n; v++) edges.push_back({stress::rand_int(0, v - 1), v, stress::rand_int(0, max_w)});
    for (int i = 0; i < extra; i++){
        edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(0, max_w)});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    for (auto& e : edges){
        if (stress::rand_int(0, 1)) swap(e[0], e[1]);
    }

    return edges;
}

/// Checks the center, the offset range, the tree and its diameter, and returns the diameter
long long check(int n, const vector<array<long long, 3>>& edges, const vector<vector<long long>>& d){
    AbsoluteCenter g(n);
    for (int i = 0; i < (int)edges.size(); i++) assert(g.add_edge(edges[i][0], edges[i][1], edges[i][2]) == i);
    AbsoluteCenterResult r = g.solve();

    if (r.edge == -1){
        assert(edges.empty() && r.diameter == 0 && r.offset == 0 && r.tree.empty());
        return 0;
    }
    assert(0 <= r.edge && r.edge < (int)edges.size());
    assert(0 <= r.offset && r.offset <= 2 * edges[r.edge][2]);
    assert(farthest(d, edges[r.edge], r.offset) == r.diameter);
    assert(tree_diameter(n, edges, r.tree) == r.diameter);
    return r.diameter;
}

int main(){
    /// Small weights, every half-integer offset tried, and the diameter against every spanning tree
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 6);
        vector<array<long long, 3>> edges = random_graph(n, stress::rand_int(0, 4), it % 2 ? 1 : 6);
        vector<vector<long long>> d = floyd(n, edges);
        long long best = check(n, edges, d);
        assert(best == brute_center(d, edges, true));
        assert(best == brute_mdst(n, edges));
    }

    /// Larger weights and denser graphs against the tent crossings
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 12);
        vector<array<long long, 3>> edges = random_graph(n, stress::rand_int(0, 3 * n), 1000000000);
        vector<vector<long long>> d = floyd(n, edges);
        assert(check(n, edges, d) == brute_center(d, edges, false));
    }

    /// Weights near the header's limit, (n - 1) * max w just below 2^61
    for (long long it = 0; it < stress::scaled(5); it++){
        int n = stress::rand_int(30, 60);
        vector<array<long long, 3>> edges = random_graph(n, 2 * n, ((1LL << 61) - 1) / (n - 1));
        vector<vector<long long>> d = floyd(n, edges);
        assert(check(n, edges, d) == brute_center(d, edges, false));
    }

    /// Paths of maximum weight edges sit at the limit, catching an infinity too small for doubled distances
    for (int n = 2; n <= 6; n++){
        vector<array<long long, 3>> limit;
        for (int v = 1; v < n; v++) limit.push_back({v - 1, v, ((1LL << 61) - 1) / (n - 1)});
        vector<vector<long long>> d = floyd(n, limit);
        assert(check(n, limit, d) == brute_center(d, limit, false));
    }

    /// Larger random graphs, and a long path whose shortest path tree is n deep
    vector<array<long long, 3>> edges = random_graph(300, 300, 1000000000), path;
    check(300, edges, floyd(300, edges));
    for (int v = 1; v < 400; v++) path.push_back({v - 1, v, stress::rand_int(0, 1000000000)});
    check(400, path, floyd(400, path));

    stress::expect_assert("graphs/absolute_center.cpp", []{
        AbsoluteCenter g(3);
        g.add_edge(0, 1, 5);
        g.solve();
    });
    stress::expect_assert("graphs/absolute_center.cpp", []{ AbsoluteCenter(2).add_edge(0, 1, -1); });

    return 0;
}
