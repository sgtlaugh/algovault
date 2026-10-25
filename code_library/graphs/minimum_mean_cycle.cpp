/***
 *
 * Minimum Mean Cycle (Karp's algorithm)
 * A directed cycle minimizing (total weight) / (number of edges), its mean exact as a reduced fraction
 *
 * Complexity: O(n (n + m)) time, O(n^2) memory
 *
 * MinimumMeanCycle g(n): empty graph on vertices 0 to n - 1
 * g.add_edge(u, v, w): directed edge u -> v, any long long weight, self loops and parallel edges allowed
 * g.solve(): {exists, num, den, cycle}, the minimum mean num / den as a reduced fraction with den > 0
 *
 * Notes:
 *   - cycle lists a simple cycle with that mean in traversal order, as edge ids numbered from 0 in add_edge order
 *   - exists is false when the graph has no cycle, then num = 0, den = 1 and cycle is empty
 *   - n * max |w| must stay <= 4e18
 *   - For the maximum mean cycle, negate the weights and the resulting num
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MeanCycle{
    bool exists;
    long long num, den;
    vector<int> cycle;
};

struct MinimumMeanCycle{
    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;

    MinimumMeanCycle(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    MeanCycle solve() const{
        const long long INF = numeric_limits<long long>::max();
        int m = edges.size();
        vector<vector<long long>> d(n + 1, vector<long long>(n, INF));
        vector<vector<int>> par(n + 1, vector<int>(n, -1));

        fill(d[0].begin(), d[0].end(), 0);
        for (int k = 1; k <= n; k++){
            for (int i = 0; i < m; i++){
                auto [u, v, w] = edges[i];
                if (d[k - 1][u] != INF && d[k - 1][u] + w < d[k][v]) d[k][v] = d[k - 1][u] + w, par[k][v] = i;
            }
        }

        int best = -1;
        long long num = 0, den = 1;
        for (int v = 0; v < n; v++){
            if (d[n][v] == INF) continue;

            /// Suffixes of the n-edge walk ending at v are walks ending at v, so every d[k][v] is finite here
            long long p = d[n][v] - d[0][v], q = n;
            for (int k = 1; k < n; k++){
                long long a = d[n][v] - d[k][v], b = n - k;
                if ((__int128)a * q > (__int128)p * b) p = a, q = b;
            }
            if (best == -1 || (__int128)p * den < (__int128)num * q) best = v, num = p, den = q;
        }
        if (best == -1) return {false, 0, 1, {}};

        long long g = gcd(num, den);
        num /= g, den /= g;

        /// Walking back n edges visits n + 1 vertices, so a vertex repeats before par[0] is ever read
        vector<int> pos(n, -1), ids(n + 1), cycle;
        int x = best, k = n;
        while (pos[x] == -1){
            pos[x] = k, ids[k] = par[k][x];
            x = edges[ids[k]].u, k--;
        }
        for (int i = k + 1; i <= pos[x]; i++) cycle.push_back(ids[i]);
        return {true, num, den, cycle};
    }
};

int main(){
    /// 0 -> 1 -> 0 has mean 3 / 2 and beats the triangle 0 -> 1 -> 2 -> 0 with mean 6 / 3
    MinimumMeanCycle g(3);
    g.add_edge(0, 1, 1);
    g.add_edge(1, 0, 2);
    g.add_edge(1, 2, 2);
    g.add_edge(2, 0, 3);
    MeanCycle res = g.solve();
    assert(res.exists);
    assert(res.num == 3);
    assert(res.den == 2);
    sort(res.cycle.begin(), res.cycle.end());
    assert((res.cycle == vector<int>{0, 1}));  /// edge ids, not vertices

    g.add_edge(2, 2, -1);
    res = g.solve();
    assert(res.num == -1);
    assert(res.den == 1);
    assert((res.cycle == vector<int>{4}));  /// a self loop is a cycle of one edge

    MinimumMeanCycle dag(3);
    dag.add_edge(0, 1, -5);
    dag.add_edge(1, 2, -5);
    dag.add_edge(0, 2, -5);
    assert(!dag.solve().exists);
    return 0;
}
