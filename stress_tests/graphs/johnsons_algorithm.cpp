#include "../common.h"

#define main library_main
#include "../../code_library/graphs/johnsons_algorithm.cpp"
#undef main

/// Floyd Warshall with explicit unreachable markers, negative cycle iff some vertex reaches itself below zero
bool floyd(int n, const vector<array<long long, 3>>& edges, vector<vector<long long>>& d){
    const long long INF = Johnson::INF;
    d.assign(n, vector<long long>(n, INF));
    for (int v = 0; v < n; v++) d[v][v] = 0;
    for (auto [u, v, w] : edges) d[u][v] = min(d[u][v], w);

    for (int k = 0; k < n; k++){
        for (int i = 0; i < n; i++){
            if (d[i][k] == INF) continue;
            for (int j = 0; j < n; j++){
                if (d[k][j] != INF) d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
            }
        }
    }

    for (int v = 0; v < n; v++){
        if (d[v][v] < 0) return false;
    }
    return true;
}

/// Solves after every batch on one instance, so stale state from an earlier solve() would show up
int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 15), m = stress::rand_int(0, 4 * n);
        bool allow_cycles = it % 4 == 0;
        vector<long long> potential(n);
        for (auto& x : potential) x = stress::rand_int(-1000000, 1000000);

        Johnson g(n);
        vector<array<long long, 3>> edges;
        int batches = stress::rand_int(1, 3);
        for (int b = 0; b < batches; b++){
            for (int i = 0; i < m; i++){
                int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
                long long w = allow_cycles ? stress::rand_int(-10, 30) : stress::rand_int(0, 1000) + potential[v] - potential[u];
                edges.push_back({u, v, w}), g.add_edge(u, v, w);
            }

            vector<vector<long long>> expected;
            bool ok = g.solve(), floyd_ok = floyd(n, edges, expected);
            assert(ok == floyd_ok);
            if (ok) assert(g.dist == expected);
            else assert(g.dist.empty());
        }
    }

    /// A long chain of negative edges, Bellman Ford needs all n rounds here
    int n = 400;
    Johnson chain(n);
    for (int i = 0; i + 1 < n; i++) chain.add_edge(i + 1, i, -1000000);
    assert(chain.solve());
    assert(chain.dist[n - 1][0] == -1000000LL * (n - 1) && chain.dist[0][n - 1] == Johnson::INF);

    /// A heavy negative self loop drops h by 1e17 per round, without the early exit it overflows long before round n
    Johnson heavy(200);
    heavy.add_edge(0, 0, -100000000000000000LL);
    assert(!heavy.solve() && heavy.dist.empty());

    /// Reweighted distances reach twice the true ones: 1 -> 2 at -B pulls h[2] down to -B, so 0 -> 2 at B becomes 2B
    const long long B = 1e18 - 1;
    Johnson wide(3);
    wide.add_edge(0, 2, B), wide.add_edge(1, 2, -B);
    assert(wide.solve() && wide.dist[0][2] == B && wide.dist[1][2] == -B);

    Johnson empty(0);
    assert(empty.solve() && empty.dist.empty());
    return 0;
}
