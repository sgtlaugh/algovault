/***
 *
 * Minimum Mean Cycle (Karp's algorithm)
 * Finds a directed cycle minimizing (total weight) / (number of edges), exact as a reduced fraction, and the cycle itself
 *
 * Complexity: O(n (n + m)) time, O(n^2) memory for the walk table
 *
 * MinimumMeanCycle g(n): n vertices, 0-based
 * g.add_edge(u, v, w): directed edge u -> v, long long weight (negative, self loops and parallel edges allowed)
 * g.solve(): returns false if the graph has no cycle, otherwise sets
 *     num / den: the minimum mean, reduced, den > 0
 *     cycle: edge ids (in add_edge order) of a simple cycle with that mean, in traversal order
 * n * max|w| must stay <= 4e18
 *
 * d[k][v] is the lightest walk of exactly k edges ending at v (starting anywhere), Karp's theorem gives
 * mean = min over v of max over k < n of (d[n][v] - d[k][v]) / (n - k)
 * The n-edge walk to the minimizing v is a shortest walk under weights w - mean, so every cycle on it has mean exactly mean
 * For maximum mean cycle, negate the weights and the resulting num
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MinimumMeanCycle{
    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;
    long long num = 0, den = 1;
    vector<int> cycle;

    MinimumMeanCycle(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    bool solve(){
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
        cycle.clear();
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
        if (best == -1) return false;

        long long g = gcd(num, den);
        num /= g, den /= g;

        /// Walking back n edges visits n + 1 vertices, so a vertex repeats before par[0] is ever read
        vector<int> pos(n, -1), ids(n + 1);
        int x = best, k = n;
        while (pos[x] == -1){
            pos[x] = k, ids[k] = par[k][x];
            x = edges[ids[k]].u, k--;
        }
        for (int i = k + 1; i <= pos[x]; i++) cycle.push_back(ids[i]);
        return true;
    }
};

int main(){
    MinimumMeanCycle g(3);
    g.add_edge(0, 1, 1);
    g.add_edge(1, 0, 2);
    g.add_edge(1, 2, 2);
    g.add_edge(2, 0, 3);
    assert(g.solve());
    assert(g.num == 3 && g.den == 2);  /// 0 -> 1 -> 0 has mean 3 / 2, the triangle has mean 6 / 3
    vector<int> two = g.cycle;
    sort(two.begin(), two.end());
    assert((two == vector<int>{0, 1}));  /// edge ids, not vertices

    g.add_edge(2, 2, -1);
    assert(g.solve());
    assert(g.num == -1 && g.den == 1);
    assert((g.cycle == vector<int>{4}));  /// a self loop is a cycle of one edge

    MinimumMeanCycle dag(3);
    dag.add_edge(0, 1, -5);
    dag.add_edge(1, 2, -5);
    dag.add_edge(0, 2, -5);
    assert(!dag.solve());
    return 0;
}
