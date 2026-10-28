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
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<pair<long long, int>>> heap;
        dist[s] = 0, heap.push({0, s});
        while (!heap.empty()){
            long long d = heap.top().first;
            int u = heap.top().second;
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
            long long total_flow = points.back().first, total_cost = points.back().second;
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
constexpr long long MCMF::INF;

int main(){
    const int dx[] = {0, 0, -1, 1}, dy[] = {-1, 1, 0, 0};
    int len;

    while (scanf("%d", &len) == 1 && len){
        vector<string> wall(5);
        static char row[1024];
        for (auto& s : wall) scanf("%s", row), s = row;

        auto in = [&](int i, int j){ return 1 + 2 * (j * 5 + i); };
        int src = 0, sink = 10 * len + 1;
        MCMF g(10 * len + 2);
        for (int j = 0; j < len; j++){
            for (int i = 0; i < 5; i++){
                int c = in(i, j);
                g.add_edge(c, c + 1, 1, wall[i][j] == '@' ? 0 : wall[i][j] - '0');
                for (int d = 0; d < 4; d++){
                    int x = i + dx[d], y = j + dy[d];
                    if (x >= 0 && x < 5 && y >= 0 && y < len) g.add_edge(c + 1, in(x, y), 1, 0);
                }
                if (j == 0 && wall[i][0] == '@') g.add_edge(src, c, 1, 0);
                if (j == len - 1) g.add_edge(c + 1, sink, 1, 0);
            }
        }
        printf("%lld\n", g.solve(src, sink).second);
    }

    return 0;
}
