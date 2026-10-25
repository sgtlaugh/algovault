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

const int MAX_BUDGET = 10000;

int main(){
    int t;

    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, r;
        scanf("%d %d", &n, &r);
        FloydWarshall fw(n);
        char row[128];
        for (int i = 0; i < n; i++){
            scanf("%s", row);
            for (int j = 0; j < n; j++){
                if (row[j] == 'Y') fw.add_edge(i, j, 1);
            }
        }
        fw.solve();

        vector<vector<double>> afford(n, vector<double>(MAX_BUDGET + 1, 1.0));
        for (int k = 1; k < n; k++){
            double window = 0;
            for (int m = 0; m <= MAX_BUDGET; m++){
                if (m >= 1) window += afford[k - 1][m - 1];
                if (m > r) window -= afford[k - 1][m - r - 1];
                afford[k][m] = window / r;
            }
        }

        int q;
        scanf("%d", &q);
        printf("Case %d\n", tc);
        while (q--){
            int a, b, m;
            scanf("%d %d %d", &a, &b, &m);
            long long d = fw.dist[a - 1][b - 1];
            printf("%.6f\n", d == FloydWarshall::INF ? 0.0 : afford[d][m] + 1e-10);
        }
        puts("");
    }
    return 0;
}
