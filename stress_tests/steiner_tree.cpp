#include "common.h"

#define main library_main
#include "../code_library/steiner_tree.cpp"
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

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 9), m = stress::rand_int(0, 3 * n);
        vector<array<long long, 3>> edges;
        for (int i = 0; i < m; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            if (u != v) edges.push_back({u, v, stress::rand_int(0, it % 2 ? 5 : 1000000000)});
        }
        int k = stress::rand_int(0, min(n, 5));
        vector<int> terminals(n);
        iota(terminals.begin(), terminals.end(), 0);
        shuffle(terminals.begin(), terminals.end(), stress::rng());
        terminals.resize(k);
        assert(steiner_tree(n, edges, terminals) == brute(n, edges, terminals));
    }

    /// Larger graph, many terminals: the answer is at most an MST and at least the farthest terminal pair distance
    int n = 3000;
    vector<array<long long, 3>> edges;
    for (int i = 1; i < n; i++) edges.push_back({i, stress::rand_int(0, i - 1), stress::rand_int(1, 1000)});
    for (int i = 0; i < 6000; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(1, 1000)});
    vector<int> terminals = {0, 7, 100, 999, 1500, 2999, 42, 1234};
    long long tree = steiner_tree(n, edges, terminals);
    long long all = 0;
    for (int i = 0; i + 1 < n; i++) all += edges[i][2];
    assert(0 < tree && tree <= all);
    assert(tree == steiner_tree(n, edges, vector<int>(terminals.rbegin(), terminals.rend())));
    return 0;
}
