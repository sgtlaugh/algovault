#include "../common.h"

#define main library_main
#include "../../code_library/graphs/graph_matching.cpp"
#undef main

/// Maximum matching by always deciding the lowest unmatched node: skip it or pair it with a neighbour
int brute(int n, const vector<int>& adj_mask){
    vector<int> memo(1 << n, -1);
    function<int(int)> best = [&](int used){
        if (used == (1 << n) - 1) return 0;
        int& res = memo[used];
        if (res != -1) return res;
        int x = __builtin_ctz(~used);
        res = best(used | 1 << x);
        for (int y = 0; y < n; y++){
            if ((adj_mask[x] >> y & 1) && !(used >> y & 1)) res = max(res, 1 + best(used | 1 << x | 1 << y));
        }
        return res;
    };
    return best(0);
}

/// The rank test errs with probability at most n / mod ~ 1e-8 per call and only ever underestimates, but a nightly
/// run makes ~1e5 calls, so a short answer is redrawn twice before it counts as a bug
int matching(Graph& g, int want){
    int got = g.maximum_matching();
    for (int retry = 0; retry < 2 && got < want; retry++) got = g.maximum_matching();
    return got;
}

int main(){
    for (long long it = 0; it < stress::scaled(1200); it++){
        int n = stress::rand_int(1, 14), density = stress::rand_int(0, 100);
        Graph g(n);
        g.rng.seed(stress::rng()());
        vector<int> adj_mask(n);

        /// Edges arrive one by one with a query after each, the matrix must be rebuilt from scratch every time
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) if (stress::rand_int(0, 99) < density) edges.push_back({u, v});
        }
        shuffle(edges.begin(), edges.end(), stress::rng());
        assert(g.maximum_matching() == 0);
        for (size_t i = 0; i < edges.size(); i++){
            auto [u, v] = edges[i];
            g.add_edge(u, v);
            adj_mask[u] |= 1 << v, adj_mask[v] |= 1 << u;
            if (i % 3 == 0 || i + 1 == edges.size()){
                int want = brute(n, adj_mask);
                assert(matching(g, want) == want);
            }
        }
    }
    return 0;
}
