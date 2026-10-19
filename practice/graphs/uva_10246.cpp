/***
 *
 * Floyd Warshall
 * All pairs shortest paths on a dense directed graph, negative edges and negative cycles allowed
 *
 * Complexity: O(n^3) time, O(n^2) memory
 *
 * FloydWarshall fw(n); fw.add_edge(u, v, w); fw.solve();
 *     0-based, directed, add both directions for an undirected edge (a negative undirected edge is a negative cycle)
 *     parallel edges keep the lightest, a negative self loop is a negative cycle, a non-negative one is ignored
 *     after solve(), dist[u][v] is the shortest distance, INF if v is unreachable from u,
 *     or NEG_INF if some u -> v walk passes through a negative cycle (no shortest path exists)
 *     path(u, v) returns the vertices of a shortest path u, ..., v, empty when dist[u][v] is INF or NEG_INF
 *     n * max |w| must stay below 2^61: walks that loop a negative cycle are clamped at NEG_INF, so they never overflow
 *
 * Call solve() once, after all edges are added
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FloydWarshall{
    static constexpr long long INF = 1LL << 62;
    static constexpr long long NEG_INF = -INF;

    int n;
    vector<vector<long long>> dist;
    vector<vector<int>> nxt;

    FloydWarshall(int n) : n(n), dist(n, vector<long long>(n, INF)), nxt(n, vector<int>(n, -1)){
        for (int v = 0; v < n; v++) dist[v][v] = 0, nxt[v][v] = v;
    }

    void add_edge(int u, int v, long long w){
        if (w < dist[u][v]) dist[u][v] = w, nxt[u][v] = v;
    }

    vector<int> path(int u, int v) const{
        if (dist[u][v] == INF || dist[u][v] == NEG_INF) return {};

        vector<int> res = {u};
        while (u != v){
            u = nxt[u][v];
            res.push_back(u);
        }
        return res;
    }

    void solve(){
        for (int k = 0; k < n; k++){
            for (int i = 0; i < n; i++){
                if (dist[i][k] == INF) continue;
                for (int j = 0; j < n; j++){
                    if (dist[k][j] == INF) continue;
                    long long d = max(dist[i][k] + dist[k][j], NEG_INF);
                    if (d < dist[i][j]) dist[i][j] = d, nxt[i][j] = nxt[i][k];
                }
            }
        }

        /// dist[k][k] < 0 exactly when k lies on a negative closed walk, which can be looped forever, so every walk through k can go arbitrarily low
        for (int k = 0; k < n; k++){
            if (dist[k][k] >= 0) continue;
            for (int i = 0; i < n; i++){
                if (dist[i][k] == INF) continue;
                for (int j = 0; j < n; j++){
                    if (dist[k][j] != INF) dist[i][j] = NEG_INF;
                }
            }
        }
    }
};

constexpr long long FloydWarshall::INF;
constexpr long long FloydWarshall::NEG_INF;

int main(){
    int n, m, q, tc = 0;
    while (scanf("%d %d %d", &n, &m, &q) == 3 && (n || m || q)){
        vector<int> cost(n);
        for (auto& x: cost) scanf("%d", &x);
        vector<array<int, 3>> edges(m);
        for (auto& e: edges) scanf("%d %d %d", &e[0], &e[1], &e[2]), e[0]--, e[1]--;

        /// Fixing the costliest feast city k, the route may only use cities no costlier than k: one Floyd Warshall per k
        const long long INF = FloydWarshall::INF;
        vector<vector<long long>> res(n, vector<long long>(n, INF));
        for (int k = 0; k < n; k++){
            FloydWarshall fw(n);
            for (auto e: edges){
                int u = e[0], v = e[1], w = e[2];
                if (cost[u] <= cost[k] && cost[v] <= cost[k]) fw.add_edge(u, v, w), fw.add_edge(v, u, w);
            }
            fw.solve();
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if (cost[i] <= cost[k] && cost[j] <= cost[k] && fw.dist[i][j] != INF){
                        res[i][j] = min(res[i][j], fw.dist[i][j] + cost[k]);
                    }
                }
            }
        }

        if (tc) puts("");
        printf("Case #%d\n", ++tc);
        while (q--){
            int u, v;
            scanf("%d %d", &u, &v);
            printf("%lld\n", res[u - 1][v - 1] == INF ? -1 : res[u - 1][v - 1]);
        }
    }
    return 0;
}
