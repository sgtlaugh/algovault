/***
 *
 * Johnson's Algorithm
 * All pairs shortest paths on a sparse directed graph with negative edges
 *
 * Complexity: O(n m) for Bellman Ford once, then O(n^2 + n m log n) for n runs of Dijkstra, O(n^2) memory for dist
 *
 * Johnson g(n); g.add_edge(u, v, w): directed edge, 0-based, long long weight (negative allowed)
 * g.solve(): returns false if the graph has a negative cycle (g.dist is left empty)
 *     otherwise g.dist[u][v] is the shortest distance, or Johnson::INF if v is unreachable from u
 *     more edges may be added and solve() called again
 *
 * Every |path weight| must stay below 1e18: reweighted distances reach up to twice that
 *
 * Bellman Ford from a virtual source gives potentials h, every edge is reweighted to w + h[u] - h[v] >= 0,
 * so plain Dijkstra works from each source and the true distance is recovered as d - h[u] + h[v]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Johnson{
    static constexpr long long INF = numeric_limits<long long>::max() / 4;

    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;
    vector<vector<long long>> dist;

    Johnson(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    bool solve(){
        dist.clear();
        vector<long long> h;
        if (!potentials(h)) return false;

        vector<vector<pair<int, long long>>> adj(n);
        for (auto [u, v, w] : edges) adj[u].push_back({v, w + h[u] - h[v]});

        dist.assign(n, vector<long long>(n, INF));
        for (int s = 0; s < n; s++) dijkstra(s, adj, h);
        return true;
    }

    void dijkstra(int s, const vector<vector<pair<int, long long>>>& adj, const vector<long long>& h){
        vector<long long>& d = dist[s];
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        d[s] = 0, heap.push({0, s});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != d[u]) continue;
            for (auto [v, w] : adj[u]){
                if (du + w < d[v]) d[v] = du + w, heap.push({d[v], v});
            }
        }

        for (int v = 0; v < n; v++){
            if (d[v] != INF) d[v] += h[v] - h[s];
        }
    }

    /// Bellman Ford from a virtual source that reaches every vertex with weight 0
    bool potentials(vector<long long>& h) const{
        h.assign(n, 0);
        for (int round = 0; round < n; round++){
            bool changed = false;
            for (auto [u, v, w] : edges){
                if (h[u] + w < h[v]){
                    h[v] = h[u] + w, changed = true;
                    /// No simple path goes this low, so it is a negative cycle, stopping now keeps h from overflowing in later rounds
                    if (h[v] < -INF) return false;
                }
            }
            if (!changed) break;
            if (round == n - 1) return false;  /// shortest paths need at most n - 1 passes, a change in pass n means a negative cycle
        }
        return true;
    }
};

int main(){
    Johnson a(3);
    a.add_edge(0, 1, 2), a.add_edge(1, 2, -15), a.add_edge(0, 2, -10);
    assert(a.solve());
    assert(a.dist[0][1] == 2 && a.dist[0][2] == -13 && a.dist[1][2] == -15 && a.dist[0][0] == 0);
    assert(a.dist[2][0] == Johnson::INF && a.dist[1][0] == Johnson::INF);

    a.add_edge(2, 0, 20);
    assert(a.solve() && a.dist[2][0] == 20 && a.dist[1][0] == 5 && a.dist[2][1] == 22);
    a.add_edge(2, 0, 10);
    assert(!a.solve() && a.dist.empty());

    Johnson cycle(3);
    cycle.add_edge(0, 1, 1), cycle.add_edge(1, 2, -3), cycle.add_edge(2, 0, 1);
    assert(!cycle.solve() && cycle.dist.empty());

    Johnson negative_loop(1), positive_loop(1);
    negative_loop.add_edge(0, 0, -1), positive_loop.add_edge(0, 0, 5);
    assert(!negative_loop.solve());
    assert(positive_loop.solve() && positive_loop.dist[0][0] == 0);

    Johnson empty(0);
    assert(empty.solve() && empty.dist.empty());

    Johnson b(4);
    b.add_edge(0, 1, 5), b.add_edge(0, 1, 3), b.add_edge(1, 2, -2), b.add_edge(2, 3, 4), b.add_edge(3, 1, -2);
    assert(b.solve());
    assert(b.dist[0][3] == 5 && b.dist[3][2] == -4 && b.dist[2][1] == 2);

    const long long B = 1e18 - 1;
    Johnson z(3);
    z.add_edge(1, 0, -B), z.add_edge(2, 1, 0), z.add_edge(0, 2, B);
    assert(z.solve());
    assert(z.dist[1][0] == -B && z.dist[0][2] == B && z.dist[2][0] == -B && z.dist[1][2] == 0 && z.dist[2][2] == 0);

    Johnson huge_cycle(2), tiny_cycle(2);
    huge_cycle.add_edge(0, 1, -B), huge_cycle.add_edge(1, 0, B - 1);
    tiny_cycle.add_edge(0, 1, -1), tiny_cycle.add_edge(1, 0, 0);
    assert(!huge_cycle.solve() && !tiny_cycle.solve());

    return 0;
}
