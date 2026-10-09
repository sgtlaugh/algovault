#include "../common.h"

#define main library_main
#include "../../code_library/graphs/minimum_path_cover.cpp"
#undef main

/// Maximum bipartite matching by DP over subsets of used right vertices
int brute_matching(int n, const vector<vector<int>>& adj){
    vector<int> best(1 << n, -1);
    best[0] = 0;
    for (int u = 0; u < n; u++){
        vector<int> next = best;
        for (int mask = 0; mask < (1 << n); mask++){
            if (best[mask] < 0) continue;
            for (int v : adj[u]){
                if (!(mask >> v & 1)) next[mask | 1 << v] = max(next[mask | 1 << v], best[mask] + 1);
            }
        }
        best = next;
    }
    return *max_element(best.begin(), best.end());
}

void check(int n, const vector<pair<int, int>>& edges){
    DAGPathCover g(n);
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges) g.add_edge(u, v), adj[u].push_back(v);

    vector<vector<char>> reach(n, vector<char>(n, 0));
    vector<vector<int>> closed(n);
    for (int s = 0; s < n; s++){
        vector<int> queue = {s};
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : adj[queue[i]]){
                if (!reach[s][v]) reach[s][v] = 1, queue.push_back(v), closed[s].push_back(v);
            }
        }
    }

    assert(g.min_disjoint_path_cover() == n - brute_matching(n, adj));
    int cover = g.min_path_cover();
    assert(cover == n - brute_matching(n, closed));

    int largest = 0;
    for (int mask = 0; mask < (1 << n); mask++){
        bool antichain = true;
        for (int u = 0; u < n && antichain; u++){
            for (int v = 0; v < n && antichain; v++){
                if ((mask >> u & 1) && (mask >> v & 1) && reach[u][v]) antichain = false;
            }
        }
        if (antichain) largest = max(largest, __builtin_popcount(mask));
    }

    auto chain = g.max_antichain();
    assert((int)chain.size() == largest && largest == cover);
    for (int u : chain){
        for (int v : chain) assert(!reach[u][v]);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(0, 11), m = n ? stress::rand_int(0, 2 * n) : 0;
        vector<int> label(n);
        iota(label.begin(), label.end(), 0);
        shuffle(label.begin(), label.end(), stress::rng());

        vector<pair<int, int>> edges;
        for (int i = 0; i < m; i++){
            int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 1);
            if (a == b) continue;
            if (a > b) swap(a, b);
            edges.push_back({label[a], label[b]});
        }
        check(n, edges);
    }

    /// Larger DAGs: Dilworth's equality must still hold and the antichain must be one
    for (long long it = 0; it < stress::scaled(5); it++){
        int n = 600;
        DAGPathCover g(n);
        vector<vector<int>> adj(n);
        for (int i = 0; i < 3000; i++){
            int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 1);
            if (a == b) continue;
            if (a > b) swap(a, b);
            g.add_edge(a, b), adj[a].push_back(b);
        }

        auto chain = g.max_antichain();
        assert((int)chain.size() == g.min_path_cover());
        assert(g.min_path_cover() <= g.min_disjoint_path_cover());

        vector<char> in(n, 0);
        for (int v : chain) in[v] = 1;
        for (int s : chain){
            vector<int> queue = {s};
            vector<char> seen(n, 0);
            for (int i = 0; i < (int)queue.size(); i++){
                for (int v : adj[queue[i]]){
                    assert(!in[v]);
                    if (!seen[v]) seen[v] = 1, queue.push_back(v);
                }
            }
        }
    }
    return 0;
}
