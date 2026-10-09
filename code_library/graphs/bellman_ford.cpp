/***
 *
 * Bellman Ford
 * Single source shortest paths on a directed graph with negative edges, and negative cycle retrieval
 *
 * Complexity: O(n m) per call, O(n + m) memory
 *
 * BellmanFord g(n); g.add_edge(u, v, w);   0-based, directed, long long weights (negative allowed)
 * g.shortest_paths(src)[v]:
 *     BellmanFord::INF if v is unreachable from src
 *     BellmanFord::NEG_INF if a negative cycle reachable from src reaches v (the distance is unbounded below)
 *     otherwise the shortest distance from src to v
 * g.negative_cycle(): the vertices of one negative cycle anywhere in the graph, not only one reachable from a source
 *     in edge order, cycle[i] -> cycle[i + 1] -> ... -> cycle[0], empty if the graph has no negative cycle
 *
 * Requires n * max|w| <= 9e18: a value below -(n - 1) * max|w| is lower than any simple path,
 * so the vertex is marked NEG_INF (or its parent chain holds a cycle) at once instead of being relaxed further
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct BellmanFord{
    struct Edge{
        int u, v;
        long long w;
    };

    static constexpr long long INF = numeric_limits<long long>::max();
    static constexpr long long NEG_INF = numeric_limits<long long>::min();

    int n;
    vector<Edge> edges;

    BellmanFord(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    vector<long long> shortest_paths(int src) const{
        long long low = -(n - 1) * max_abs_weight();
        vector<long long> dist(n, INF);
        dist[src] = 0;

        /// Unaffected vertices are exact after n - 1 passes, any edge still relaxing after that leads into a negative cycle's reach
        for (int pass = 0; pass < n; pass++){
            bool changed = false;
            for (auto [u, v, w] : edges){
                if (dist[u] == INF || dist[u] == NEG_INF || dist[u] + w >= dist[v]) continue;
                dist[v] = (pass == n - 1 || dist[u] + w < low) ? NEG_INF : dist[u] + w;
                changed = true;
            }
            if (!changed) break;
        }

        for (bool changed = true; changed;){
            changed = false;
            for (auto [u, v, w] : edges){
                if (dist[u] == NEG_INF && dist[v] != NEG_INF) dist[v] = NEG_INF, changed = true;
            }
        }
        return dist;
    }

    vector<int> negative_cycle() const{
        long long low = -(n - 1) * max_abs_weight();
        vector<long long> dist(n, 0);
        vector<int> parent(n, -1);

        for (int pass = 0; pass < n; pass++){
            int relaxed = -1;
            for (auto [u, v, w] : edges){
                if (dist[u] + w >= dist[v]) continue;
                dist[v] = dist[u] + w, parent[v] = u, relaxed = v;
                /// Below every simple path weight, so the parent chain of v cannot end at a root: it holds a cycle
                if (dist[v] < low) return cycle_through(parent, v);
            }
            if (relaxed == -1) return {};
            /// An acyclic parent chain of at most n - 1 edges would bound dist[relaxed] by its value after n - 1 passes, so the chain holds a cycle
            if (pass == n - 1) return cycle_through(parent, relaxed);
        }
        return {};
    }

    long long max_abs_weight() const{
        long long res = 0;
        for (auto [u, v, w] : edges) res = max(res, abs(w));
        return res;
    }

    vector<int> cycle_through(const vector<int>& parent, int v) const{
        for (int i = 0; i < n; i++) v = parent[v];

        vector<int> cycle = {v};
        for (int u = parent[v]; u != v; u = parent[u]) cycle.push_back(u);
        reverse(cycle.begin(), cycle.end());
        return cycle;
    }
};

bool is_rotation(vector<int> cycle, const vector<int>& expected){
    for (size_t i = 0; i < cycle.size(); i++){
        if (cycle == expected) return true;
        rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
    }
    return false;
}

int main(){
    const long long INF = BellmanFord::INF, NEG_INF = BellmanFord::NEG_INF;

    BellmanFord g(6);
    g.add_edge(0, 1, 4), g.add_edge(0, 2, 2), g.add_edge(2, 1, -1), g.add_edge(1, 3, 3), g.add_edge(3, 4, -2);
    assert((g.shortest_paths(0) == vector<long long>{0, 1, 2, 4, 2, INF}));
    assert((g.shortest_paths(3) == vector<long long>{INF, INF, INF, 0, -2, INF}));
    assert(g.negative_cycle().empty());

    BellmanFord h(5);
    h.add_edge(0, 1, 1), h.add_edge(1, 2, -1), h.add_edge(2, 1, -1), h.add_edge(2, 3, 5), h.add_edge(4, 0, 1);
    assert((h.shortest_paths(0) == vector<long long>{0, NEG_INF, NEG_INF, NEG_INF, INF}));
    assert((h.shortest_paths(4) == vector<long long>{1, NEG_INF, NEG_INF, NEG_INF, 0}));
    assert((h.shortest_paths(3) == vector<long long>{INF, INF, INF, 0, INF}));
    assert(is_rotation(h.negative_cycle(), {1, 2}));

    BellmanFord zero(2);
    zero.add_edge(0, 1, -3), zero.add_edge(1, 0, 3);
    assert((zero.shortest_paths(0) == vector<long long>{0, -3}));
    assert(zero.negative_cycle().empty());

    BellmanFord loop(1);
    loop.add_edge(0, 0, -1);
    assert((loop.shortest_paths(0) == vector<long long>{NEG_INF}));
    assert((loop.negative_cycle() == vector<int>{0}));

    BellmanFord hidden(3);
    hidden.add_edge(1, 2, -5), hidden.add_edge(2, 1, 2);
    assert((hidden.shortest_paths(0) == vector<long long>{0, INF, INF}));
    assert(is_rotation(hidden.negative_cycle(), {1, 2}));

    BellmanFord triangle(3);
    triangle.add_edge(2, 0, 1), triangle.add_edge(1, 2, -1), triangle.add_edge(0, 1, -1);
    assert(is_rotation(triangle.negative_cycle(), {0, 1, 2}));
    assert((triangle.shortest_paths(1) == vector<long long>{NEG_INF, NEG_INF, NEG_INF}));

    const long long W = 3000000000000000000LL;
    BellmanFord heavy(3);
    heavy.add_edge(0, 1, -W), heavy.add_edge(1, 2, -W);
    assert((heavy.shortest_paths(0) == vector<long long>{0, -W, -2 * W}));
    assert(heavy.negative_cycle().empty());
    heavy.add_edge(2, 0, -W);
    assert((heavy.shortest_paths(0) == vector<long long>{NEG_INF, NEG_INF, NEG_INF}));
    assert(is_rotation(heavy.negative_cycle(), {0, 1, 2}));

    assert(BellmanFord(0).negative_cycle().empty());

    return 0;
}
