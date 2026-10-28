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
    const double SCALE = 1e9;
    const long long EPS = 1000;
    int n, m;

    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        vector<int> x(n), y(n);
        for (int i = 0; i < n; i++) scanf("%d %d", &x[i], &y[i]);
        auto dist2 = [&](int i, int j){ return (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]); };
        auto length = [&](int i, int j){ return llround(sqrt((double)dist2(i, j)) * SCALE); };

        FloydWarshall fw(n);
        vector<vector<bool>> road(n, vector<bool>(n, false));
        for (int i = 0; i < m; i++){
            int a, b;
            scanf("%d %d", &a, &b);
            a--, b--;
            road[a][b] = road[b][a] = true;
            fw.add_edge(a, b, length(a, b)), fw.add_edge(b, a, length(a, b));
        }
        fw.solve();
        const vector<vector<long long>>& d = fw.dist;

        int bu = -1, bv = -1;
        long long best = SCALE;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (road[u][v]) continue;
                long long w = length(u, v), gain = 0;
                for (int i = 0; i < n; i++){
                    for (int j = 0; j < n; j++) gain += d[i][j] - min({d[i][j], d[i][u] + w + d[v][j], d[i][v] + w + d[u][j]});
                }

                if (gain > best + EPS) best = gain, bu = u, bv = v;
                else if (bu != -1 && gain >= best - EPS && dist2(u, v) < dist2(bu, bv)) bu = u, bv = v;
            }
        }

        if (bu == -1) puts("No road required");
        else printf("%d %d\n", bu + 1, bv + 1);
    }

    return 0;
}
