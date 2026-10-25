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
            for (Edge e : edges){
                int u = e.u, v = e.v;
                long long w = e.w;
                if (dist[u] == INF || dist[u] == NEG_INF || dist[u] + w >= dist[v]) continue;
                dist[v] = (pass == n - 1 || dist[u] + w < low) ? NEG_INF : dist[u] + w;
                changed = true;
            }
            if (!changed) break;
        }

        for (bool changed = true; changed;){
            changed = false;
            for (Edge e : edges){
                int u = e.u, v = e.v;
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
            for (Edge e : edges){
                int u = e.u, v = e.v;
                long long w = e.w;
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
        for (Edge e : edges) res = max(res, abs(e.w));
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

constexpr long long BellmanFord::INF;

bool is_rotation(vector<int> cycle, const vector<int>& expected){
    for (size_t i = 0; i < cycle.size(); i++){
        if (cycle == expected) return true;
        rotate(cycle.begin(), cycle.begin() + 1, cycle.end());
    }
    return false;
}

int main(){
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, m;
        char a[32], b[32];
        scanf("%d", &n);
        map<string, int> id;
        for (int i = 0; i < n; i++) scanf("%s", a), id[a] = i;

        BellmanFord g(n * n);
        scanf("%d", &m);
        while (m--){
            int w;
            scanf("%s %s %d", a, b, &w);
            int u = id[a], v = id[b];
            for (int h = 0; h + 1 < n; h++) g.add_edge(h * n + u, (h + 1) * n + v, w);
        }
        vector<long long> dist = g.shortest_paths(0);

        if (tc > 1) puts("");
        printf("Scenario #%d\n", tc);
        int q;
        scanf("%d", &q);
        while (q--){
            int stops;
            scanf("%d", &stops);
            long long res = BellmanFord::INF;
            for (int h = 1; h <= min(stops + 1, n - 1); h++) res = min(res, dist[h * n + n - 1]);
            if (res == BellmanFord::INF) puts("No satisfactory flights");
            else printf("Total cost of flight(s) is $%lld\n", res);
        }
    }
    return 0;
}
