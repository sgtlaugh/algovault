#include "../common.h"

#define main library_main
#include "../../code_library/graphs/steiner_tree.cpp"
#undef main

/// Minimum over vertex subsets containing every terminal of the MST of the induced subgraph, when connected
long long brute(int n, const vector<array<long long, 3>>& edges, const vector<int>& terminals){
    if (terminals.size() <= 1) return 0;
    int need = 0;
    for (int t : terminals) need |= 1 << t;

    vector<int> order(edges.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){ return edges[a][2] < edges[b][2]; });

    long long best = STEINER_INF;
    for (int mask = 0; mask < (1 << n); mask++){
        if ((mask & need) != need) continue;
        vector<int> parent(n);
        iota(parent.begin(), parent.end(), 0);
        function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };
        long long total = 0;
        int joined = 0;
        for (int i : order){
            auto [u, v, w] = edges[i];
            if (!(mask >> u & 1) || !(mask >> v & 1)) continue;
            int a = find(u), b = find(v);
            if (a != b) parent[a] = b, total += w, joined++;
        }
        if (joined == __builtin_popcount(mask) - 1) best = min(best, total);
    }
    return best;
}

/// Same subset recurrence, but spreading through Floyd-Warshall distances instead of Dijkstra, for graphs too big for brute()
/// It shares the recurrence with the library, so brute() on n <= 9 is the independent oracle; this one guards the Dijkstra spreading
long long floyd_reference(int n, const vector<array<long long, 3>>& edges, const vector<int>& terminals){
    int k = terminals.size();
    if (k <= 1) return 0;

    vector<vector<long long>> dist(n, vector<long long>(n, STEINER_INF));
    for (int v = 0; v < n; v++) dist[v][v] = 0;
    for (auto [u, v, w] : edges){
        dist[u][v] = min(dist[u][v], w);
        dist[v][u] = min(dist[v][u], w);
    }
    for (int m = 0; m < n; m++){
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++) dist[u][v] = min(dist[u][v], dist[u][m] + dist[m][v]);
        }
    }

    vector<vector<long long>> dp(1 << k, vector<long long>(n, STEINER_INF));
    for (int i = 0; i < k; i++) dp[1 << i] = dist[terminals[i]];
    for (int mask = 1; mask < (1 << k); mask++){
        if (__builtin_popcount(mask) == 1) continue;
        vector<long long> merged(n, STEINER_INF);
        for (int sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask){
            for (int v = 0; v < n; v++) merged[v] = min(merged[v], dp[sub][v] + dp[mask ^ sub][v]);
        }
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++) dp[mask][v] = min(dp[mask][v], merged[u] + dist[u][v]);
        }
    }

    return *min_element(dp[(1 << k) - 1].begin(), dp[(1 << k) - 1].end());
}

/// Tiny weights force ties, 1e15 weights catch any narrowing of the 64-bit weight type
/// Bit 31 stays clear so a narrowed weight is wrong but nonnegative: an assertion fires instead of Dijkstra looping
long long random_weight(long long it){
    if (it % 3 == 2) return stress::rand_int(0, 1000000000000000LL) & ~(1LL << 31);
    return stress::rand_int(0, it % 3 ? 1000000000 : 5);
}

vector<int> random_terminals(int n, int max_k){
    vector<int> terminals(n);
    iota(terminals.begin(), terminals.end(), 0);
    shuffle(terminals.begin(), terminals.end(), stress::rng());
    terminals.resize(stress::rand_int(0, min(n, max_k)));
    return terminals;
}

/// Small graphs against the subset brute force, medium ones against the Floyd-Warshall recurrence,
/// several terminal sets per graph to check solve() leaves the graph intact
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 9), m = stress::rand_int(0, 3 * n);
        vector<array<long long, 3>> edges;
        SteinerTree g(n);
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = random_weight(it);
            edges.push_back({u, v, w});
            g.add_edge(u, v, w);
        }

        for (int q = 0; q < 3; q++){
            vector<int> terminals = random_terminals(n, 5);
            assert(g.solve(terminals) == brute(n, edges, terminals));
        }
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 40), m = stress::rand_int(0, 3 * n);
        vector<array<long long, 3>> edges;
        SteinerTree g(n);
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = random_weight(it);
            edges.push_back({u, v, w});
            g.add_edge(u, v, w);
        }

        vector<int> terminals = random_terminals(n, 6);
        assert(g.solve(terminals) == floyd_reference(n, edges, terminals));
    }

    /// Larger graph, many terminals: the answer is positive, at most the spanning tree weight, and independent of terminal order
    int n = 3000;
    SteinerTree g(n);
    long long all = 0;
    for (int i = 1; i < n; i++){
        long long w = stress::rand_int(1, 1000);
        g.add_edge(i, stress::rand_int(0, i - 1), w);
        all += w;
    }
    for (int i = 0; i < 6000; i++) g.add_edge(stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(1, 1000));
    vector<int> terminals = {0, 7, 100, 999, 1500, 2999, 42, 1234};

    long long tree = g.solve(terminals);
    assert(0 < tree && tree <= all);
    assert(tree == g.solve(vector<int>(terminals.rbegin(), terminals.rend())));

    return 0;
}
