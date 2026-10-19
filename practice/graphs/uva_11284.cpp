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

int cents(){
    double x;
    scanf("%lf", &x);
    return (int)lround(x * 100);
}

int main(){
    int t, n, m, p;
    scanf("%d", &t);
    while (t--){
        scanf("%d %d", &n, &m);
        FloydWarshall fw(n + 1);
        for (int i = 0; i < m; i++){
            int a, b;
            scanf("%d %d", &a, &b);
            int c = cents();
            fw.add_edge(a, b, c), fw.add_edge(b, a, c);
        }
        fw.solve();

        /// Operas at the same store are bought together, so only the distinct stores enter the bitmask
        map<int, long long> saving;
        scanf("%d", &p);
        for (int i = 0; i < p; i++){
            int a;
            scanf("%d", &a);
            saving[a] += cents();
        }
        vector<int> stores;
        vector<long long> gain;
        for (auto p : saving) stores.push_back(p.first), gain.push_back(p.second);

        int k = stores.size();
        vector<vector<long long>> dp(1 << k, vector<long long>(k, FloydWarshall::INF));
        for (int i = 0; i < k; i++) dp[1 << i][i] = fw.dist[0][stores[i]];
        long long best = 0;
        for (int mask = 1; mask < (1 << k); mask++){
            long long total = 0;
            for (int i = 0; i < k; i++){
                if (mask >> i & 1) total += gain[i];
            }
            for (int i = 0; i < k; i++){
                if (dp[mask][i] == FloydWarshall::INF) continue;
                best = max(best, total - dp[mask][i] - fw.dist[stores[i]][0]);
                for (int j = 0; j < k; j++){
                    if (mask >> j & 1) continue;
                    long long& cell = dp[mask | 1 << j][j];
                    cell = min(cell, dp[mask][i] + fw.dist[stores[i]][stores[j]]);
                }
            }
        }

        if (best > 0) printf("Daniel can save $%lld.%02lld\n", best / 100, best % 100);
        else puts("Don't leave the house");
    }
    return 0;
}
