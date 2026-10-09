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

/// Blossom's pairs must be distinct edges of the graph, as many as the claimed size
void check_mates(const Blossom& b, const vector<int>& adj_mask, int size){
    int pairs = 0;
    for (int v = 0; v < b.n; v++){
        int u = b.mate[v];
        if (u == -1) continue;
        assert(u >= 0 && u < b.n && u != v && b.mate[u] == v && (adj_mask[v] >> u & 1));
        pairs += u > v;
    }
    assert(pairs == size);
}

/// Same check on graphs too large for a bitmask, against a sorted edge list
void check_mates(const Blossom& b, const vector<pair<int, int>>& sorted_edges, int size){
    int pairs = 0;
    for (int v = 0; v < b.n; v++){
        int u = b.mate[v];
        if (u == -1) continue;
        assert(u >= 0 && u < b.n && u != v && b.mate[u] == v);
        assert(binary_search(sorted_edges.begin(), sorted_edges.end(), make_pair(min(u, v), max(u, v))));
        pairs += u > v;
    }
    assert(pairs == size);
}

/// Every graph on up to 6 vertices, Blossom against the brute force
void exhaustive_small(){
    for (int n = 1; n <= 6; n++){
        vector<pair<int, int>> slots;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) slots.push_back({u, v});
        }

        for (int mask = 0; mask < 1 << slots.size(); mask++){
            Blossom b(n);
            vector<int> adj_mask(n);
            for (size_t i = 0; i < slots.size(); i++){
                if (!(mask >> i & 1)) continue;
                auto [u, v] = slots[i];
                b.add_edge(u, v);
                adj_mask[u] |= 1 << v, adj_mask[v] |= 1 << u;
            }

            int want = brute(n, adj_mask), got = b.maximum_matching();
            assert(got == want);
            check_mates(b, adj_mask, got);
        }
    }
}

/// Larger graphs, Blossom (given some edges twice) against the Tutte rank: random sparse ones, odd cycles chained by bridges (nested blossoms),
/// and up to 400 vertices with about one edge per vertex, where augmenting paths run long through big blossoms
void random_large(){
    for (long long it = 0; it < stress::scaled(450); it++){
        int n = it % 3 == 2 ? stress::rand_int(100, 400) : stress::rand_int(2, 120);
        vector<pair<int, int>> edges;
        if (it % 3 == 1){
            for (int start = 0; start < n;){
                int len = min(n - start, 2 * (int)stress::rand_int(1, 4) + 1);
                for (int i = 0; i < len; i++){
                    int u = start + i, v = start + (i + 1) % len;
                    if (u != v) edges.push_back({u, v});
                }
                if (start > 0) edges.push_back({stress::rand_int(0, start - 1), start + stress::rand_int(0, len - 1)});
                start += len;
            }
        }
        else{
            int m = it % 3 == 0 ? stress::rand_int(0, 3 * n) : stress::rand_int(n, 3 * n / 2);
            for (int i = 0; i < m; i++){
                int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
                if (u != v) edges.push_back({u, v});
            }
        }

        /// Random labels, so the greedy initial pairs do not follow the generator's structure
        vector<int> label(n);
        iota(label.begin(), label.end(), 0);
        shuffle(label.begin(), label.end(), stress::rng());
        for (auto& [u, v] : edges){
            u = label[u], v = label[v];
            if (u > v) swap(u, v);
        }

        sort(edges.begin(), edges.end());
        edges.erase(unique(edges.begin(), edges.end()), edges.end());
        auto order = edges;
        shuffle(order.begin(), order.end(), stress::rng());
        Graph g(n);
        g.rng.seed(stress::rng()());
        Blossom b(n);
        for (auto [u, v] : order){
            g.add_edge(u, v);
            if (stress::rand_int(0, 1)) swap(u, v);
            b.add_edge(u, v);
            if (stress::rand_int(0, 9) == 0) b.add_edge(v, u);
        }

        int got = b.maximum_matching();
        assert(matching(g, got) == got);
        check_mates(b, edges, got);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(1200); it++){
        int n = stress::rand_int(1, 14), density = stress::rand_int(0, 100);
        Graph g(n);
        g.rng.seed(stress::rng()());
        Blossom b(n);
        vector<int> adj_mask(n);

        /// Edges arrive one by one with a query after each, both structs must rebuild from scratch every time
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) if (stress::rand_int(0, 99) < density) edges.push_back({u, v});
        }
        shuffle(edges.begin(), edges.end(), stress::rng());

        assert(g.maximum_matching() == 0 && b.maximum_matching() == 0);
        for (size_t i = 0; i < edges.size(); i++){
            auto [u, v] = edges[i];
            g.add_edge(u, v);
            b.add_edge(u, v);
            adj_mask[u] |= 1 << v, adj_mask[v] |= 1 << u;
            if (i % 3 == 0 || i + 1 == edges.size()){
                int want = brute(n, adj_mask);
                assert(matching(g, want) == want);
                assert(b.maximum_matching() == want);
                check_mates(b, adj_mask, want);
            }
        }
    }

    exhaustive_small();
    random_large();
    return 0;
}
