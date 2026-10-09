#include "../common.h"

#define main library_main
#include "../../code_library/graphs/bellman_ford.cpp"
#undef main

using Edges = vector<array<long long, 3>>;

const __int128 FLOYD_INF = (__int128)1 << 120;

/// Floyd Warshall in __int128: v is NEG_INF from s iff s reaches some k with d[k][k] < 0 that reaches v
vector<vector<__int128>> floyd(int n, const Edges& edges){
    vector<vector<__int128>> d(n, vector<__int128>(n, FLOYD_INF));
    for (int v = 0; v < n; v++) d[v][v] = 0;
    for (auto [u, v, w] : edges) d[u][v] = min(d[u][v], (__int128)w);

    for (int k = 0; k < n; k++){
        for (int i = 0; i < n; i++){
            if (d[i][k] == FLOYD_INF) continue;
            for (int j = 0; j < n; j++){
                if (d[k][j] != FLOYD_INF) d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
            }
        }
    }
    return d;
}

void check(int n, const Edges& edges){
    BellmanFord g(n);
    for (auto [u, v, w] : edges) g.add_edge(u, v, w);
    auto d = floyd(n, edges);

    for (int s = 0; s < n; s++){
        vector<long long> got = g.shortest_paths(s);
        for (int v = 0; v < n; v++){
            bool unbounded = false;
            for (int k = 0; k < n; k++) unbounded |= d[s][k] != FLOYD_INF && d[k][k] < 0 && d[k][v] != FLOYD_INF;
            if (unbounded) assert(got[v] == BellmanFord::NEG_INF);
            else if (d[s][v] == FLOYD_INF) assert(got[v] == BellmanFord::INF);
            else assert(got[v] == d[s][v]);
        }
    }

    bool has_cycle = false;
    for (int k = 0; k < n; k++) has_cycle |= d[k][k] < 0;
    vector<int> cycle = g.negative_cycle();
    assert(has_cycle == !cycle.empty());
    if (!has_cycle) return;

    vector<vector<__int128>> best(n, vector<__int128>(n, FLOYD_INF));
    for (auto [u, v, w] : edges) best[u][v] = min(best[u][v], (__int128)w);
    set<int> distinct(cycle.begin(), cycle.end());
    assert(distinct.size() == cycle.size());
    __int128 total = 0;
    for (size_t i = 0; i < cycle.size(); i++){
        int u = cycle[i], v = cycle[(i + 1) % cycle.size()];
        assert(0 <= u && u < n && best[u][v] != FLOYD_INF);
        total += best[u][v];
    }
    assert(total < 0);
}

int main(){
    /// Every graph on up to 3 vertices with each ordered pair (self loops too) absent, -1 or +1, so zero cycles show up as well
    for (int n = 1; n <= 3; n++){
        int pairs = n * n, total = 1;
        for (int i = 0; i < pairs; i++) total *= 3;
        for (int mask = 0; mask < total; mask++){
            Edges edges;
            for (int i = 0, rest = mask; i < pairs; i++, rest /= 3){
                if (rest % 3) edges.push_back({i / n, i % n, rest % 3 == 1 ? -1 : 1});
            }
            shuffle(edges.begin(), edges.end(), stress::rng());
            check(n, edges);
        }
    }

    /// Random multigraphs: weights at the n * max|w| = 9e18 boundary, potential-shifted weights (negative edges, no negative cycle), small weights
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 9), m = stress::rand_int(0, 3 * n);
        long long w_max = 9000000000000000000LL / n;
        vector<long long> potential(n);
        for (auto& p : potential) p = stress::rand_int(-1000, 1000);

        Edges edges;
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = stress::rand_int(-30, 30);
            if (it % 3 == 0) w = stress::rand_int(0, 1) ? stress::rand_int(-w_max, w_max) : (stress::rand_int(0, 1) ? w_max : -w_max);
            if (it % 3 == 1) w = stress::rand_int(0, 50) + potential[v] - potential[u];
            edges.push_back({u, v, w});
        }
        check(n, edges);
    }

    /// A chain listed backwards needs all n - 1 passes, then a heavy negative cycle at its end poisons everything after it
    const int n = 1500;
    const long long w = 9000000000000000000LL / n;
    BellmanFord chain(n);
    for (int i = n - 2; i >= 0; i--) chain.add_edge(i, i + 1, -w);
    vector<long long> got = chain.shortest_paths(0);
    for (int v = 0; v < n; v++) assert(got[v] == -w * v);
    assert(chain.negative_cycle().empty());

    int start = n / 2;
    chain.add_edge(start + 10, start, -w);
    got = chain.shortest_paths(0);
    for (int v = 0; v < n; v++) assert(got[v] == (v < start ? -w * v : BellmanFord::NEG_INF));
    vector<int> cycle = chain.negative_cycle();
    vector<int> expected(11);
    iota(expected.begin(), expected.end(), start);
    assert(is_rotation(cycle, expected));

    return 0;
}
