#include "common.h"

#define main library_main
#include "../code_library/johnsons_algorithm.cpp"
#undef main

/// Floyd Warshall with explicit unreachable markers, negative cycle iff some vertex reaches itself below zero
bool floyd(int n, const vector<array<long long, 3>>& edges, vector<vector<long long>>& d){
    const long long INF = JOHNSON_INF;
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

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 15), m = stress::rand_int(0, 4 * n);
        bool allow_cycles = it % 4 == 0;
        vector<long long> potential(n);
        for (auto& x : potential) x = stress::rand_int(-1000000, 1000000);
        vector<array<long long, 3>> edges;
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = allow_cycles ? stress::rand_int(-10, 30) : stress::rand_int(0, 1000) + potential[v] - potential[u];
            edges.push_back({u, v, w});
        }

        vector<vector<long long>> got, expected;
        bool ok = johnson(n, edges, got), floyd_ok = floyd(n, edges, expected);
        assert(ok == floyd_ok);
        if (ok) assert(got == expected);
        else assert(got.empty());
    }

    /// A long chain of negative edges, Bellman Ford needs all n rounds here
    int n = 400;
    vector<array<long long, 3>> chain;
    for (int i = 0; i + 1 < n; i++) chain.push_back({i + 1, i, -1000000});
    vector<vector<long long>> got;
    assert(johnson(n, chain, got));
    assert(got[n - 1][0] == -1000000LL * (n - 1) && got[0][n - 1] == JOHNSON_INF);

    /// A heavy negative self loop drops h by 1e17 per round, without the early exit it overflows long before round n
    assert(!johnson(200, {{0, 0, -100000000000000000LL}}, got) && got.empty());
    return 0;
}
