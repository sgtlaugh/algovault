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

/// The library's subset recurrence spread by Floyd-Warshall instead of Dijkstra: it guards only the spreading
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

/// Minimum weight edge subset joining each pair, over all 2^m subsets: independent of the vertex subset brute force
long long brute_forest(int n, const vector<array<long long, 3>>& edges, const vector<int>& terminals){
    long long best = STEINER_INF;
    for (int mask = 0; mask < (1 << edges.size()); mask++){
        vector<int> parent(n);
        iota(parent.begin(), parent.end(), 0);
        function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };
        long long total = 0;
        for (int i = 0; i < (int)edges.size(); i++){
            if (mask >> i & 1) parent[find(edges[i][0])] = find(edges[i][1]), total += edges[i][2];
        }

        bool joined = true;
        for (int i = 0; i + 1 < (int)terminals.size(); i += 2) joined &= find(terminals[i]) == find(terminals[i + 1]);
        if (joined) best = min(best, total);
    }

    return best;
}

/// Asserts the ids are distinct edges joining every terminal, and returns their total weight
long long tree_weight(int n, const vector<array<long long, 3>>& edges, const vector<int>& terminals, vector<int> ids){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };
    vector<bool> seen(edges.size(), false);
    long long total = 0;
    for (int id : ids){
        assert(0 <= id && id < (int)edges.size() && !seen[id]);
        seen[id] = true;
        parent[find(edges[id][0])] = find(edges[id][1]);
        total += edges[id][2];
    }

    for (int t : terminals) assert(find(t) == find(terminals[0]));
    return total;
}

/// Tiny weights force ties, 1e15 weights catch any narrowing of the 64-bit weight type
long long random_weight(long long it){
    /// a narrowed weight stays nonnegative, so an assertion fires instead of Dijkstra looping
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

/// Several terminal sets per graph check that solve() leaves the graph intact
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
            int k = terminals.size();
            vector<long long> best = g.solve_subsets(terminals);
            assert((int)best.size() == 1 << k);
            for (int mask = 0; mask < (1 << k); mask++){
                vector<int> subset;
                for (int i = 0; i < k; i++){
                    if (mask >> i & 1) subset.push_back(terminals[i]);
                }
                assert(best[mask] == brute(n, edges, subset));
            }
            assert(g.solve(terminals) == best.back());
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
        long long best = floyd_reference(n, edges, terminals);
        assert(g.solve(terminals) == best);
        if (best < STEINER_INF) assert(tree_weight(n, edges, terminals, g.tree(terminals)) == best);
    }

    /// tree against every edge subset, each terminal paired with the first so brute_forest asks for one component
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 7), m = stress::rand_int(0, 10);
        vector<array<long long, 3>> edges;
        SteinerTree g(n);
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = random_weight(it);
            edges.push_back({u, v, w});
            assert(g.add_edge(u, v, w) == i);
        }

        vector<int> terminals = random_terminals(n, 5), pairs;
        for (int t : terminals) pairs.push_back(terminals[0]), pairs.push_back(t);
        long long best = brute_forest(n, edges, pairs);
        if (best == STEINER_INF) continue;
        vector<int> ids = g.tree(terminals);
        assert(tree_weight(n, edges, terminals, ids) == best);
        if (terminals.size() <= 1) assert(ids.empty());
    }

    stress::expect_assert("graphs/steiner_tree.cpp", []{ SteinerTree(3).tree({0, 2}); });

    /// Steiner forest by the header's recipe, terminals paired as (0, 1), (2, 3), (4, 5)
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(2, 7), m = stress::rand_int(0, 11);
        vector<array<long long, 3>> edges;
        SteinerTree g(n);
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long w = random_weight(it);
            edges.push_back({u, v, w});
            g.add_edge(u, v, w);
        }

        vector<int> terminals = random_terminals(n, 6);
        terminals.resize(terminals.size() / 2 * 2);
        int pairs = terminals.size() / 2;
        vector<long long> best = g.solve_subsets(terminals);
        vector<long long> forest(1 << pairs, STEINER_INF);
        forest[0] = 0;
        for (int s = 1; s < (1 << pairs); s++){
            for (int t = s; t > 0; t = (t - 1) & s){
                int covered = 0;
                for (int i = 0; i < pairs; i++){
                    if (t >> i & 1) covered |= 3 << (2 * i);
                }
                forest[s] = min(forest[s], forest[s ^ t] + best[covered]);
            }
        }
        assert(forest.back() == brute_forest(n, edges, terminals));
    }

    /// Larger graph, many terminals: positive, at most the spanning tree weight, independent of order, tree agrees
    int n = 3000;
    vector<array<long long, 3>> edges;
    SteinerTree g(n);
    long long all = 0;
    for (int i = 0; i < n - 1 + 6000; i++){
        int u = i < n - 1 ? i + 1 : stress::rand_int(0, n - 1);
        int v = i < n - 1 ? stress::rand_int(0, i) : stress::rand_int(0, n - 1);
        long long w = stress::rand_int(1, 1000);
        edges.push_back({u, v, w});
        g.add_edge(u, v, w);
        if (i < n - 1) all += w;
    }
    vector<int> terminals = {0, 7, 100, 999, 1500, 2999, 42, 1234};

    long long tree = g.solve(terminals);
    assert(0 < tree && tree <= all);
    assert(tree == g.solve(vector<int>(terminals.rbegin(), terminals.rend())));
    assert(tree_weight(n, edges, terminals, g.tree(terminals)) == tree);

    return 0;
}
