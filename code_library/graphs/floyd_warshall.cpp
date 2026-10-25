/***
 *
 * Floyd Warshall
 * All pairs shortest paths on a dense directed graph, negative edges and negative cycles allowed
 *
 * Complexity:
 *   - constructor: O(n^2) time and memory
 *   - relax_through: O(n^2)
 *   - solve: O(n^3)
 *
 * FloydWarshall fw(n): vertices 0 to n - 1, no edges
 * fw.add_edge(u, v, w): directed edge, parallel edges keep the lightest
 * fw.relax_through(k): every dist[u][v] may now use k as an inner vertex, assumes no negative cycle
 * fw.solve(): computes every distance, call it once after all edges are added
 * fw.path(u, v): the vertices u, ..., v of a shortest path, empty when dist[u][v] is INF or NEG_INF
 * fw.dist[u][v]: the shortest distance, INF if unreachable, NEG_INF if some u -> v walk passes a negative cycle
 *
 * Notes:
 *   - Add both directions for an undirected edge, so a negative undirected edge is a negative cycle
 *   - A negative self loop is a negative cycle, a non-negative one is ignored
 *   - n * max |w| must stay below 2^61
 *   - After relax_through(k) for each k in a set S, dist[u][v] is the shortest u -> v path with all inner vertices in S
 *   - So adding vertices one at a time in any order gives the distances among the vertices added so far
 *   - After solve(), add_edge(u, v, w) then relax_through(u) and relax_through(v) updates every dist in O(n^2)
 *   - path stays correct after relax_through
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

    void relax_through(int k){
        for (int i = 0; i < n; i++){
            if (dist[i][k] == INF) continue;
            for (int j = 0; j < n; j++){
                if (dist[k][j] == INF) continue;
                long long d = max(dist[i][k] + dist[k][j], NEG_INF);
                if (d < dist[i][j]) dist[i][j] = d, nxt[i][j] = nxt[i][k];
            }
        }
    }

    void solve(){
        for (int k = 0; k < n; k++) relax_through(k);

        /// dist[k][k] < 0 exactly when k lies on a negative closed walk, so every walk through k goes arbitrarily low
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

    vector<int> path(int u, int v) const{
        if (dist[u][v] == INF || dist[u][v] == NEG_INF) return {};

        vector<int> res = {u};
        while (u != v){
            u = nxt[u][v];
            res.push_back(u);
        }

        return res;
    }
};

int main(){
    const long long INF = FloydWarshall::INF, NEG_INF = FloydWarshall::NEG_INF;

    /// 0 -> 2 -> 1 beats the direct 4, and nothing leaves 3
    FloydWarshall a(4);
    a.add_edge(0, 1, 4), a.add_edge(0, 2, 1), a.add_edge(2, 1, 2), a.add_edge(1, 3, 1), a.add_edge(2, 3, 5);
    a.solve();
    assert(a.dist[0][1] == 3);
    assert(a.dist[3][0] == INF);
    assert((a.path(0, 3) == vector<int>{0, 2, 1, 3}));

    /// 1 <-> 2 is a negative cycle: everything reachable through it has no shortest path
    FloydWarshall b(4);
    b.add_edge(0, 1, 1), b.add_edge(1, 2, -1), b.add_edge(2, 1, -1), b.add_edge(2, 3, 2);
    b.solve();
    assert(b.dist[0][3] == NEG_INF);
    assert(b.dist[3][3] == 0);  /// 3 is past the cycle, its own distance is unaffected
    assert(b.path(0, 3).empty());

    /// Vertices join in the order 2, 0, 1, and each relax_through lets paths pass through one more of them
    FloydWarshall c(3);
    c.add_edge(0, 1, 5), c.add_edge(0, 2, 1), c.add_edge(2, 1, 2), c.add_edge(1, 0, 1);
    c.relax_through(2);
    assert(c.dist[0][1] == 3);
    assert(c.dist[1][2] == INF);
    c.relax_through(0);
    assert(c.dist[1][2] == 2);
    c.relax_through(1);
    assert(c.dist[2][0] == 3);

    /// A new edge 3 -> 0 of weight 1 after solve: 3 now reaches everything through 0
    a.add_edge(3, 0, 1);
    a.relax_through(3), a.relax_through(0);
    assert(a.dist[3][1] == 4);
    assert((a.path(3, 1) == vector<int>{3, 0, 2, 1}));
    assert(a.dist[2][0] == 4);
    return 0;
}
