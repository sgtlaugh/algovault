/***
 *
 * Min Cost Max Flow with Dijkstra and potentials
 * Successive shortest paths where Johnson potentials keep reduced costs non-negative
 *
 * Complexity: O(F m log n) for total flow F, plus O(m) per call, O(n m) if a residual edge has negative cost
 *
 * MCMF g(n): empty graph on nodes 0 to n - 1
 * g.add_edge(u, v, cap, cost): adds a directed edge and returns its index, add both directions for undirected edges
 * g.slope(s, t, limit): augments like solve, returns the breakpoints {flow, cost} of min cost as a function of flow
 * g.solve(s, t, limit): pushes up to limit units (default: as much as possible), returns {flow, cost}
 * g.flow(id): flow on edge id after solve or slope
 *
 * Notes:
 *   - slope() starts at {0, 0}, ends at what solve returns, linear between breakpoints, slopes strictly increasing
 *   - Repeated solve or slope calls continue from the residual graph, results count only the newly added flow and cost
 *   - Costs may be negative, but a negative cycle reachable from s makes solve and slope loop forever
 *   - flow * |cost| summed over a path must fit in long long
 *   - The default min cost flow, mcmf_spfa.cpp has the same API minus slope() and is faster for small integer costs
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MCMF{
    static constexpr long long INF = numeric_limits<long long>::max() / 4;

    struct Edge{
        int to;
        long long cap, cost;
    };

    int n;
    vector<Edge> edges;
    vector<vector<int>> adj;
    vector<long long> potential, dist;
    vector<int> parent_edge;

    MCMF(int n) : n(n), adj(n), potential(n, 0), dist(n), parent_edge(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        adj[u].push_back(edges.size()), edges.push_back({v, cap, cost});
        adj[v].push_back(edges.size()), edges.push_back({u, 0, -cost});
        return edges.size() - 2;
    }

    vector<pair<long long, long long>> slope(int s, int t, long long limit = INF){
        assert(s != t);
        init_potential(s);
        vector<pair<long long, long long>> points = {{0, 0}};
        long long prev_unit_cost = 0;

        while (points.back().first < limit && dijkstra(s, t)){
            auto [total_flow, total_cost] = points.back();
            long long push = limit - total_flow, unit_cost = 0;
            for (int v = t; v != s; v = edges[parent_edge[v] ^ 1].to){
                push = min(push, edges[parent_edge[v]].cap);
                unit_cost += edges[parent_edge[v]].cost;
            }
            for (int v = t; v != s; v = edges[parent_edge[v] ^ 1].to){
                edges[parent_edge[v]].cap -= push;
                edges[parent_edge[v] ^ 1].cap += push;
            }

            if (points.size() > 1 && unit_cost == prev_unit_cost) points.pop_back();
            points.push_back({total_flow + push, total_cost + push * unit_cost});
            prev_unit_cost = unit_cost;
        }
        return points;
    }

    pair<long long, long long> solve(int s, int t, long long limit = INF){
        return slope(s, t, limit).back();
    }

    long long flow(int id) const{
        return edges[id ^ 1].cap;
    }

    /// Zero potentials work without negative residual costs; else SPFA, as index-order Bellman Ford can need n passes
    void init_potential(int s){
        bool negative = false;
        for (const Edge& e : edges) negative |= e.cap > 0 && e.cost < 0;
        fill(potential.begin(), potential.end(), negative ? INF : 0);
        if (!negative) return;

        vector<char> queued(n, 0);
        queue<int> q;
        potential[s] = 0, q.push(s), queued[s] = 1;
        while (!q.empty()){
            int u = q.front();
            q.pop(), queued[u] = 0;
            for (int id : adj[u]){
                const Edge& e = edges[id];
                if (e.cap > 0 && potential[u] + e.cost < potential[e.to]){
                    potential[e.to] = potential[u] + e.cost;
                    if (!queued[e.to]) q.push(e.to), queued[e.to] = 1;
                }
            }
        }
    }

    bool dijkstra(int s, int t){
        fill(dist.begin(), dist.end(), INF);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        dist[s] = 0, heap.push({0, s});
        while (!heap.empty()){
            auto [d, u] = heap.top();
            heap.pop();
            if (d != dist[u]) continue;
            for (int id : adj[u]){
                const Edge& e = edges[id];
                if (e.cap <= 0) continue;
                long long nd = d + e.cost + potential[u] - potential[e.to];
                if (nd < dist[e.to]){
                    dist[e.to] = nd, parent_edge[e.to] = id;
                    heap.push({nd, e.to});
                }
            }
        }

        /// Only reachable nodes get new potentials, unreachable ones can never become reachable again
        for (int u = 0; u < n; u++){
            if (dist[u] < INF) potential[u] += dist[u];
        }
        return dist[t] < INF;
    }
};

int main(){
    /***
     * 0 -> 1 (cap 2, cost 1), 0 -> 2 (cap 1, cost 2), 1 -> 2 (cap 1, cost 1), 1 -> 3 (cap 1, cost 3), 2 -> 3 (cap 2, cost 1)
     * Max flow 3: paths 0-1-2-3 (cost 3), 0-2-3 (cost 3), 0-1-3 (cost 4), total cost 10
    ***/
    MCMF g(4);
    g.add_edge(0, 1, 2, 1), g.add_edge(0, 2, 1, 2), g.add_edge(1, 2, 1, 1), g.add_edge(1, 3, 1, 3);
    int e23 = g.add_edge(2, 3, 2, 1);

    MCMF limited = g, curve = g;

    assert((g.solve(0, 3) == make_pair(3LL, 10LL)));
    assert(g.flow(e23) == 2);

    assert((limited.solve(0, 3, 1) == make_pair(1LL, 3LL)));  /// one unit along a cost 3 path
    assert((limited.solve(0, 3) == make_pair(2LL, 7LL)));     /// continues: the other 2 units cost 3 + 4

    using Points = vector<pair<long long, long long>>;
    assert((curve.slope(0, 3) == Points{{0, 0}, {2, 6}, {3, 10}}));  /// both units at cost 3 share a segment

    MCMF negative(3);
    negative.add_edge(0, 1, 5, -4), negative.add_edge(1, 2, 3, 2), negative.add_edge(0, 2, 4, 1);
    assert((negative.solve(0, 2) == make_pair(7LL, -2LL)));  /// 3 units at -2 each, then 4 at 1
    return 0;
}
