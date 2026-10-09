#include "../common.h"

#define main library_main
#include "../../code_library/graphs/floyd_warshall.cpp"
#undef main

const long long INF = FloydWarshall::INF, NEG_INF = FloydWarshall::NEG_INF;

/// Bellman Ford from s: after n - 1 rounds an edge that still relaxes ends on a vertex that a negative cycle reaches,
/// and everything reachable from such a vertex is NEG_INF
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

void check(int n, const vector<array<long long, 3>>& edges){
    FloydWarshall fw(n);
    vector<vector<long long>> lightest(n, vector<long long>(n, INF));
    for (auto [u, v, w] : edges){
        fw.add_edge(u, v, w);
        lightest[u][v] = min(lightest[u][v], w);
    }
    fw.solve();

    for (int s = 0; s < n; s++){
        assert(fw.dist[s] == bellman_ford(n, edges, s));
        for (int t = 0; t < n; t++){
            vector<int> p = fw.path(s, t);
            if (fw.dist[s][t] == INF || fw.dist[s][t] == NEG_INF){
                assert(p.empty());
                continue;
            }

            assert(!p.empty() && (int)p.size() <= n && p.front() == s && p.back() == t);
            long long total = 0;
            for (int i = 0; i + 1 < (int)p.size(); i++){
                assert(lightest[p[i]][p[i + 1]] != INF);
                total += lightest[p[i]][p[i + 1]];
            }
            assert(total == fw.dist[s][t]);
        }
    }
}

int main(){
    /// Every graph on up to 3 vertices with edge weights from {absent, -1, 0, 2} and self loops from {absent, -1},
    /// a non-negative self loop never changes a distance
    const long long choices[4] = {INF, -1, 0, 2};
    for (int n = 0; n <= 3; n++){
        int cells = n * n, total = 1;
        for (int i = 0; i < cells; i++) total *= 4;
        for (int mask = 0; mask < total; mask++){
            vector<array<long long, 3>> edges;
            bool skip = false;
            for (int i = 0, x = mask; i < cells; i++, x /= 4){
                if (i / n == i % n && x % 4 > 1) skip = true;
                if (x % 4) edges.push_back({i / n, i % n, choices[x % 4]});
            }
            if (!skip) check(n, edges);
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 2500 ? stress::rand_int(1, 12) : stress::rand_int(13, 40);
        int m = stress::rand_int(0, it % 3 == 0 ? min(n * n, 150) : 3 * n);
        bool allow_cycles = it % 2 == 0;
        vector<long long> potential(n);
        for (auto& x : potential) x = stress::rand_int(-1000000000000LL, 1000000000000LL);

        vector<array<long long, 3>> edges;
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = allow_cycles ? stress::rand_int(-3, 12) : stress::rand_int(0, 2) + potential[v] - potential[u];
            edges.push_back({u, v, w});
        }
        check(n, edges);
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
