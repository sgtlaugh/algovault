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
    const long long INF = FloydWarshall::INF, NEG_INF = FloydWarshall::NEG_INF;

    FloydWarshall a(4);
    a.add_edge(0, 1, 4), a.add_edge(0, 2, 1), a.add_edge(2, 1, 2), a.add_edge(1, 3, 1), a.add_edge(2, 3, 5);
    a.solve();
    assert(a.dist[0][1] == 3 && a.dist[0][3] == 4 && a.dist[2][3] == 3 && a.dist[3][3] == 0);
    assert(a.dist[3][0] == INF && a.dist[1][2] == INF);
    assert((a.path(0, 3) == vector<int>{0, 2, 1, 3}));
    assert((a.path(2, 2) == vector<int>{2}));
    assert(a.path(3, 0).empty());

    FloydWarshall b(5);
    b.add_edge(0, 1, 1), b.add_edge(1, 2, -1), b.add_edge(2, 1, -1), b.add_edge(2, 3, 2), b.add_edge(4, 0, 3);
    b.solve();
    assert(b.dist[0][1] == NEG_INF && b.dist[0][3] == NEG_INF && b.dist[4][3] == NEG_INF && b.dist[1][1] == NEG_INF);
    assert(b.dist[4][0] == 3 && b.dist[3][3] == 0 && b.dist[3][0] == INF && b.dist[0][4] == INF && b.dist[1][0] == INF);
    assert((b.path(4, 0) == vector<int>{4, 0}));
    assert(b.path(0, 3).empty() && b.path(3, 1).empty());

    FloydWarshall c(2);
    c.add_edge(0, 0, -1), c.add_edge(0, 1, 5);
    c.solve();
    assert(c.dist[0][0] == NEG_INF && c.dist[0][1] == NEG_INF && c.dist[1][1] == 0 && c.dist[1][0] == INF);

    FloydWarshall d(1);
    d.add_edge(0, 0, 5);
    d.solve();
    assert(d.dist[0][0] == 0 && (d.path(0, 0) == vector<int>{0}));

    FloydWarshall e(3);
    e.add_edge(0, 1, 5), e.add_edge(0, 1, -2), e.add_edge(1, 2, -3), e.add_edge(0, 2, -4);
    e.solve();
    assert(e.dist[0][1] == -2 && e.dist[0][2] == -5 && e.dist[1][2] == -3);
    assert((e.path(0, 2) == vector<int>{0, 1, 2}));

    FloydWarshall f(3);
    f.add_edge(0, 1, 0), f.add_edge(1, 0, 0), f.add_edge(1, 2, 0);
    f.solve();
    assert(f.dist[0][2] == 0 && f.dist[1][0] == 0 && f.dist[2][0] == INF);
    assert((f.path(0, 2) == vector<int>{0, 1, 2}));

    FloydWarshall empty(0);
    empty.solve();
    assert(empty.dist.empty());
    return 0;
}
