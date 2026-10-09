/***
 *
 * Min Cost Max Flow with Dijkstra and potentials
 * Successive shortest paths where Johnson potentials keep reduced costs non-negative
 *
 * Complexity: O(F * m log n) for total flow F, plus O(n * m) once for the initial Bellman Ford
 *
 * MCMF g(n); g.add_edge(u, v, cap, cost): directed edge, returns its index, add both directions for undirected edges
 * g.solve(s, t, limit): pushes up to limit units (default: as much as possible), returns {flow, cost}
 * g.flow(id): flow on edge id after solve
 *
 * Costs may be negative as long as there is no negative cycle
 * flow * |cost| summed over a path must fit in long long
 *
 * The default min cost flow choice, see mcmf_spfa.cpp for an SPFA variant without potentials
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
        assert(s != t);
        init_potential(s);
        long long total_flow = 0, total_cost = 0;

        while (total_flow < limit && dijkstra(s, t)){
            long long push = limit - total_flow;
            for (int v = t; v != s; v = edges[parent_edge[v] ^ 1].to) push = min(push, edges[parent_edge[v]].cap);
            for (int v = t; v != s; v = edges[parent_edge[v] ^ 1].to){
                edges[parent_edge[v]].cap -= push;
                edges[parent_edge[v] ^ 1].cap += push;
                total_cost += push * edges[parent_edge[v]].cost;
            }
            total_flow += push;
        }
        return {total_flow, total_cost};
    }
};

int main(){
    /***
     * 0 -> 1 (cap 2, cost 1), 0 -> 2 (cap 1, cost 2), 1 -> 2 (cap 1, cost 1), 1 -> 3 (cap 1, cost 3), 2 -> 3 (cap 2, cost 1)
     * Max flow 3: paths 0-1-2-3 (cost 3), 0-2-3 (cost 3), 0-1-3 (cost 4), total cost 10
    ***/
    MCMF g(4);
    int e01 = g.add_edge(0, 1, 2, 1);
    g.add_edge(0, 2, 1, 2), g.add_edge(1, 2, 1, 1), g.add_edge(1, 3, 1, 3);
    int e23 = g.add_edge(2, 3, 2, 1);
    assert((g.solve(0, 3) == make_pair(3LL, 10LL)));
    assert(g.flow(e01) == 2 && g.flow(e23) == 2);

    MCMF limited(4);
    limited.add_edge(0, 1, 2, 1), limited.add_edge(0, 2, 1, 2), limited.add_edge(1, 2, 1, 1), limited.add_edge(1, 3, 1, 3), limited.add_edge(2, 3, 2, 1);
    assert((limited.solve(0, 3, 1) == make_pair(1LL, 3LL)));

    MCMF negative(3);
    negative.add_edge(0, 1, 5, -4), negative.add_edge(1, 2, 3, 2), negative.add_edge(0, 2, 4, 1);
    assert((negative.solve(0, 2) == make_pair(7LL, -6LL + 4LL)));

    MCMF disconnected(3);
    disconnected.add_edge(0, 1, 5, 1);
    assert((disconnected.solve(0, 2) == make_pair(0LL, 0LL)));

    return 0;
}
