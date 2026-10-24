#include "../common.h"

#define main library_main
#include "../../code_library/graphs/floyd_warshall.cpp"
#undef main

const long long INF = FloydWarshall::INF, NEG_INF = FloydWarshall::NEG_INF;

/// Bellman Ford from s: everything reachable from an edge that still relaxes after n - 1 rounds is NEG_INF
vector<long long> bellman_ford(int n, const vector<array<long long, 3>>& edges, int s){
    vector<long long> d(n, INF);
    d[s] = 0;
    for (int round = 0; round + 1 < n; round++){
        for (auto [u, v, w] : edges){
            if (d[u] != INF && d[u] + w < d[v]) d[v] = d[u] + w;
        }
    }

    vector<bool> neg(n, false);
    for (int round = 0; round < n; round++){
        for (auto [u, v, w] : edges){
            if (d[u] != INF && (neg[u] || d[u] + w < d[v])) neg[v] = true;
        }
    }

    for (int v = 0; v < n; v++){
        if (neg[v]) d[v] = NEG_INF;
    }
    return d;
}

vector<vector<long long>> lightest_edges(int n, const vector<array<long long, 3>>& edges){
    vector<vector<long long>> lightest(n, vector<long long>(n, INF));
    for (auto [u, v, w] : edges) lightest[u][v] = min(lightest[u][v], w);
    return lightest;
}

/// inner marks the vertices a path may pass through, so the walk must stay inside them apart from its ends
void check_path(const FloydWarshall& fw, const vector<vector<long long>>& lightest, const vector<bool>& inner, int s, int t){
    vector<int> p = fw.path(s, t);
    if (fw.dist[s][t] == INF || fw.dist[s][t] == NEG_INF){
        assert(p.empty());
        return;
    }

    assert(!p.empty() && (int)p.size() <= fw.n && p.front() == s && p.back() == t);
    long long total = 0;
    for (int i = 0; i + 1 < (int)p.size(); i++){
        assert(lightest[p[i]][p[i + 1]] != INF);
        total += lightest[p[i]][p[i + 1]];
        if (i > 0) assert(inner[p[i]]);
    }
    assert(total == fw.dist[s][t]);
}

void check(int n, const vector<array<long long, 3>>& edges){
    FloydWarshall fw(n);
    for (auto [u, v, w] : edges) fw.add_edge(u, v, w);
    fw.solve();

    vector<vector<long long>> lightest = lightest_edges(n, edges);
    vector<bool> inner(n, true);
    for (int s = 0; s < n; s++){
        assert(fw.dist[s] == bellman_ford(n, edges, s));
        for (int t = 0; t < n; t++) check_path(fw, lightest, inner, s, t);
    }
}

/// Reverse vertex deletion: relax_through in a random order, each prefix against a Bellman Ford limited to that prefix
void check_vertex_order(int n, const vector<array<long long, 3>>& edges){
    FloydWarshall fw(n);
    for (auto [u, v, w] : edges) fw.add_edge(u, v, w);

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    for (int i = n - 1; i > 0; i--) swap(order[i], order[stress::rand_int(0, i)]);

    vector<vector<long long>> lightest = lightest_edges(n, edges);
    vector<bool> inner(n, false);
    for (int k : order){
        fw.relax_through(k);
        inner[k] = true;
        for (int s = 0; s < n; s++){
            vector<array<long long, 3>> allowed;
            for (auto e : edges){
                if (e[0] == s || inner[e[0]]) allowed.push_back(e);
            }
            assert(fw.dist[s] == bellman_ford(n, allowed, s));
            for (int t = 0; t < n; t++) check_path(fw, lightest, inner, s, t);
        }
    }
}

/// Edges arrive one at a time after solve(), each followed by relax_through on its two ends
void check_new_edges(int n, const vector<array<long long, 3>>& edges, int initial){
    FloydWarshall fw(n);
    vector<array<long long, 3>> added(edges.begin(), edges.begin() + initial);
    for (auto [u, v, w] : added) fw.add_edge(u, v, w);
    fw.solve();

    vector<bool> inner(n, true);
    for (int e = initial; e < (int)edges.size(); e++){
        auto [u, v, w] = edges[e];
        fw.add_edge(u, v, w);
        fw.relax_through(u), fw.relax_through(v);
        added.push_back(edges[e]);

        vector<vector<long long>> lightest = lightest_edges(n, added);
        for (int s = 0; s < n; s++){
            assert(fw.dist[s] == bellman_ford(n, added, s));
            for (int t = 0; t < n; t++) check_path(fw, lightest, inner, s, t);
        }
    }
}

vector<array<long long, 3>> random_edges(int n, int m, bool allow_cycles){
    vector<long long> potential(n);
    for (auto& x : potential) x = stress::rand_int(-1000000000000LL, 1000000000000LL);

    /// Weights 0 to 2 plus a potential difference keep every cycle non-negative, zero weight cycles included
    vector<array<long long, 3>> edges;
    for (int i = 0; i < m; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        long long w = allow_cycles ? stress::rand_int(-3, 12) : stress::rand_int(0, 2) + potential[v] - potential[u];
        edges.push_back({u, v, w});
    }
    return edges;
}

int main(){
    /// Every graph on up to 3 vertices with edge weights from {absent, -1, 0, 2}
    const long long choices[4] = {INF, -1, 0, 2};
    for (int n = 0; n <= 3; n++){
        int cells = n * n, total = 1;
        for (int i = 0; i < cells; i++) total *= 4;
        for (int mask = 0; mask < total; mask++){
            vector<array<long long, 3>> edges;
            bool skip = false;
            for (int i = 0, x = mask; i < cells; i++, x /= 4){
                if (i / n == i % n && x % 4 > 1) skip = true;  /// a non-negative self loop never changes a distance
                if (x % 4) edges.push_back({i / n, i % n, choices[x % 4]});
            }
            if (!skip) check(n, edges);
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 2500 ? stress::rand_int(1, 12) : stress::rand_int(13, 40);
        int m = stress::rand_int(0, it % 3 == 0 ? min(n * n, 150) : 3 * n);
        check(n, random_edges(n, m, it % 2 == 0));
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 10), m = stress::rand_int(0, 3 * n);
        vector<array<long long, 3>> edges = random_edges(n, m, false);
        check_vertex_order(n, edges);
        check_new_edges(n, edges, stress::rand_int(0, m));
    }

    /// The header bound n * max |w| < 2^61: chains of the heaviest allowed edges in both signs
    for (int n : {2, 8, 64}){
        long long w = ((1LL << 61) - 1) / n;
        for (long long sign : {1, -1}){
            vector<array<long long, 3>> edges;
            for (int i = 0; i + 1 < n; i++) edges.push_back({i, i + 1, sign * w});
            for (int i = 0; i + 2 < n; i++) edges.push_back({i, i + 2, sign * w * 2});
            check(n, edges);
        }
    }

    /// A complete graph of heavy negative edges: unclamped walks would reach -2^60 * 2^k after k phases and overflow
    int n = 64;
    long long w = ((1LL << 61) - 1) / n;
    FloydWarshall fw(n);
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++) fw.add_edge(u, v, -w);
    }
    fw.solve();
    for (int u = 0; u < n; u++){
        for (int v = 0; v < n; v++) assert(fw.dist[u][v] == NEG_INF && fw.path(u, v).empty());
    }
    return 0;
}
