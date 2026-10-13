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

int main(){
    const int INF = 1 << 20;

    int n, m, tc = 0;
    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        FloydWarshall fw(n);
        for (int i = 0; i < m; i++){
            int a, b;
            scanf("%d %d", &a, &b);
            fw.add_edge(a - 1, b - 1, 1);
        }
        fw.solve();

        vector<vector<int>> dis(n, vector<int>(n));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) dis[i][j] = fw.dist[i][j] == FloydWarshall::INF ? INF : fw.dist[i][j];
        }

        printf("Network %d\n", ++tc);
        const int s = 0, t = 1;
        if (dis[s][t] == INF || dis[t][s] == INF){
            puts("Impossible\n");
            continue;
        }

        /// Dijkstra over pairs (i, j): moving to (k, l) pays the shortest walks j -> k -> l -> i,
        /// minus the node shared with the previous loop
        vector<vector<int>> dp(n, vector<int>(n));
        vector<vector<bool>> visited(n, vector<bool>(n, false));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) dp[i][j] = min(INF, dis[s][i] + dis[i][j] + dis[j][s] + 1);
        }

        priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>, greater<tuple<int, int, int>>> pq;
        pq.push({dp[s][s], s, s});
        while (!pq.empty()){
            auto [d, i, j] = pq.top();
            pq.pop();
            if (visited[i][j] || d != dp[i][j]) continue;
            if (i == t && j == t) break;
            visited[i][j] = true;

            for (int k = 0; k < n; k++){
                if (k == i || k == j) continue;
                for (int l = 0; l < n; l++){
                    if (l == i || l == j || visited[k][l]) continue;
                    int w = d + dis[j][k] + dis[k][l] + dis[l][i] - 1;
                    if (w < dp[k][l]) dp[k][l] = w, pq.push({w, k, l});
                }
            }
        }
        printf("Minimum number of nodes = %d\n\n", dp[t][t]);
    }
    return 0;
}
