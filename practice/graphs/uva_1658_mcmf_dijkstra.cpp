/***
 *
 * Min Cost Max Flow with Dijkstra and potentials
 * Successive shortest paths where Johnson potentials keep reduced costs non-negative
 *
 * Complexity: O(F * m log n) for total flow F, plus O(n * m) per solve or slope call for its initial Bellman Ford
 *
 * MCMF g(n); g.add_edge(u, v, cap, cost): directed edge, returns its index, add both directions for undirected edges
 * g.solve(s, t, limit): pushes up to limit units (default: as much as possible), returns {flow, cost}
 * g.slope(s, t, limit): same augmentation as solve, returns the breakpoints {flow, cost} of min cost as a function of flow
 *   starts at {0, 0}, ends at what solve returns, cost is linear between consecutive points
 *   slopes strictly increase (convex), collinear points are merged, so the min cost of any x <= max flow interpolates
 * g.flow(id): flow on edge id after solve or slope
 * Repeated solve or slope calls continue from the residual graph, results count only the newly added flow and cost
 *
 * Costs may be negative as long as there is no negative cycle
 * flow * |cost| summed over a path must fit in long long
 *
 * The default min cost flow choice, mcmf_spfa.cpp has the same API minus slope(), augments along all tied shortest paths per phase, faster when costs are small integers
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

    long long flow(int id) const{
        return edges[id ^ 1].cap;
    }

    /// Bellman Ford from s so that negative costs give valid starting potentials
    void init_potential(int s){
        fill(potential.begin(), potential.end(), INF);
        potential[s] = 0;
        for (int round = 0; round < n; round++){
            bool changed = false;
            for (int u = 0; u < n; u++){
                if (potential[u] == INF) continue;
                for (int id : adj[u]){
                    const Edge& e = edges[id];
                    if (e.cap > 0 && potential[u] + e.cost < potential[e.to]) potential[e.to] = potential[u] + e.cost, changed = true;
                }
            }
            if (!changed) break;
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

    pair<long long, long long> solve(int s, int t, long long limit = INF){
        return slope(s, t, limit).back();
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
};

int main(){
    int v, e, a, b, c;
    while (scanf("%d %d", &v, &e) == 2){
        /// Waypoint i splits into 2i -> 2i + 1 with capacity 1 so the routes share no waypoint, the ends take both ships
        MCMF g(2 * v);
        for (int i = 0; i < v; i++) g.add_edge(2 * i, 2 * i + 1, i == 0 || i == v - 1 ? 2 : 1, 0);
        while (e--) scanf("%d %d %d", &a, &b, &c), g.add_edge(2 * a - 1, 2 * b - 2, 1, c);
        printf("%lld\n", g.solve(0, 2 * v - 1, 2).second);
    }
    return 0;
}
