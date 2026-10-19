/***
 *
 * Min Cost Max Flow with SPFA
 * Successive shortest paths where each SPFA phase pushes flow along every tied shortest path with one DFS
 *
 * Complexity: exponential in theory, fast in practice, see https://codeforces.com/blog/entry/70740
 *
 * MCMFSpfa g(n); g.add_edge(u, v, cap, cost): directed edge, returns its index, add both directions for undirected edges
 * g.solve(s, t, limit): pushes up to limit units (default: as much as possible), returns {flow, cost}
 * g.flow(id): flow on edge id after solve
 * Repeated solve calls continue from the residual graph, results count only the newly added flow and cost
 * Nodes are 0 to n - 1, 1 based graphs use MCMFSpfa(n + 1)
 *
 * Costs may be negative as long as there is no negative cycle
 * flow * |cost| summed over a path must fit in long long
 * The DFS recurses once per edge of a shortest path: 1.5e5 edges overflow an 8 MB stack at -O2, use mcmf_dijkstra.cpp for long paths
 *
 * Same API as mcmf_dijkstra.cpp minus slope(), prefer this one when costs are small integers with many tied shortest paths
 *     300 x 300 assignment with costs 0 to 3: 20x faster than mcmf_dijkstra.cpp, about 2x slower on random sparse graphs
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MCMFSpfa{
    static constexpr long long INF = numeric_limits<long long>::max() / 4;

    struct Edge{
        int to;
        long long cap, cost;
    };

    int n;
    vector<Edge> edges;
    vector<vector<int>> adj;
    vector<long long> dis;
    vector<char> visited;

    MCMFSpfa(int n) : n(n), adj(n), dis(n), visited(n) {}

    int add_edge(int u, int v, long long cap, long long cost){
        adj[u].push_back(edges.size()), edges.push_back({v, cap, cost});
        adj[v].push_back(edges.size()), edges.push_back({u, 0, -cost});
        return edges.size() - 2;
    }

    long long flow(int id) const{
        return edges[id ^ 1].cap;
    }

    bool spfa(int s, int t){
        fill(dis.begin(), dis.end(), INF);
        fill(visited.begin(), visited.end(), 0);
        queue<int> q;
        dis[s] = 0, q.push(s), visited[s] = 1;

        while (!q.empty()){
            int u = q.front();
            q.pop(), visited[u] = 0;
            for (int id : adj[u]){
                const Edge& e = edges[id];
                if (e.cap > 0 && dis[u] + e.cost < dis[e.to]){
                    dis[e.to] = dis[u] + e.cost;
                    if (!visited[e.to]) q.push(e.to), visited[e.to] = 1;
                }
            }
        }

        return dis[t] < INF;
    }

    /// The sink is never marked so one DFS can push along several shortest paths
    long long dfs(int u, int t, long long lim){
        if (u == t) return lim;
        visited[u] = 1;
        long long aug = 0;

        for (int id : adj[u]){
            if (aug == lim) break;
            Edge& e = edges[id];
            if (!visited[e.to] && e.cap > 0 && dis[e.to] == dis[u] + e.cost){
                long long x = dfs(e.to, t, min(e.cap, lim - aug));
                e.cap -= x, edges[id ^ 1].cap += x, aug += x;
            }
        }

        return aug;
    }

    pair<long long, long long> solve(int s, int t, long long limit = INF){
        assert(s != t);
        long long flow = 0, cost = 0;

        while (flow < limit && spfa(s, t)){
            while (flow < limit){
                fill(visited.begin(), visited.end(), 0);
                long long aug = dfs(s, t, limit - flow);
                if (!aug) break;
                flow += aug, cost += aug * dis[t];
            }
        }

        return {flow, cost};
    }
};
constexpr long long MCMFSpfa::INF;

int main(){
    int v, e, a, b, c;
    while (scanf("%d %d", &v, &e) == 2){
        /// Waypoint i splits into 2i -> 2i + 1 with capacity 1 so the routes share no waypoint, the ends take both ships
        MCMFSpfa g(2 * v);
        for (int i = 0; i < v; i++) g.add_edge(2 * i, 2 * i + 1, i == 0 || i == v - 1 ? 2 : 1, 0);
        while (e--) scanf("%d %d %d", &a, &b, &c), g.add_edge(2 * a - 1, 2 * b - 2, 1, c);
        printf("%lld\n", g.solve(0, 2 * v - 1, 2).second);
    }
    return 0;
}
