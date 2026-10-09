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
    auto sorted_cycle = [](MinimumMeanCycle& g){
        vector<int> c = g.cycle;
        sort(c.begin(), c.end());
        return c;
    };

    MinimumMeanCycle triangle(3);
    triangle.add_edge(0, 1, 1);
    triangle.add_edge(1, 2, 2);
    triangle.add_edge(2, 0, 3);
    triangle.add_edge(1, 0, 4);
    assert(triangle.solve() && triangle.num == 2 && triangle.den == 1);
    assert((sorted_cycle(triangle) == vector<int>{0, 1, 2}));

    MinimumMeanCycle half(3);
    half.add_edge(2, 2, 2);
    half.add_edge(0, 1, 1);
    half.add_edge(1, 0, 2);
    assert(half.solve() && half.num == 3 && half.den == 2);
    assert((sorted_cycle(half) == vector<int>{1, 2}));

    MinimumMeanCycle reduced(2);
    reduced.add_edge(0, 1, 2);
    reduced.add_edge(1, 0, 4);
    assert(reduced.solve() && reduced.num == 3 && reduced.den == 1);

    MinimumMeanCycle zero(2);
    zero.add_edge(0, 1, -7);
    zero.add_edge(1, 0, 7);
    assert(zero.solve() && zero.num == 0 && zero.den == 1);

    MinimumMeanCycle negative(2);
    negative.add_edge(0, 1, -3);
    negative.add_edge(1, 0, -4);
    negative.add_edge(1, 1, -5);
    assert(negative.solve() && negative.num == -5 && negative.den == 1);
    assert((negative.cycle == vector<int>{2}));

    MinimumMeanCycle parallel(2);
    parallel.add_edge(0, 1, 10);
    parallel.add_edge(0, 1, 1);
    parallel.add_edge(1, 0, 2);
    assert(parallel.solve() && parallel.num == 3 && parallel.den == 2);
    assert((sorted_cycle(parallel) == vector<int>{1, 2}));

    MinimumMeanCycle islands(5);
    islands.add_edge(0, 1, 5);
    islands.add_edge(1, 0, 5);
    islands.add_edge(2, 3, 0);
    islands.add_edge(3, 4, 0);
    islands.add_edge(4, 2, -1);
    islands.add_edge(1, 2, -100);
    assert(islands.solve() && islands.num == -1 && islands.den == 3);
    assert((sorted_cycle(islands) == vector<int>{2, 3, 4}));

    MinimumMeanCycle dag(4);
    dag.add_edge(0, 1, -1);
    dag.add_edge(1, 2, -1);
    dag.add_edge(0, 3, -1);
    dag.add_edge(3, 2, -1);
    assert(!dag.solve() && dag.cycle.empty() && dag.num == 0 && dag.den == 1);
    dag.add_edge(2, 0, 2);
    assert(dag.solve() && dag.num == 0 && dag.den == 1 && dag.cycle.size() == 3);
    dag.add_edge(2, 2, -1);
    assert(dag.solve() && dag.num == -1 && dag.den == 1);
    assert((dag.cycle == vector<int>{5}));
    assert(!MinimumMeanCycle(0).solve());
    assert(!MinimumMeanCycle(1).solve());

    MinimumMeanCycle extreme(4);
    extreme.add_edge(0, 1, -1000000000000000000LL);
    extreme.add_edge(1, 2, -1000000000000000000LL);
    extreme.add_edge(2, 3, -1000000000000000000LL);
    extreme.add_edge(3, 0, 1000000000000000000LL);
    extreme.add_edge(3, 3, 1000000000000000000LL);
    assert(extreme.solve() && extreme.num == -500000000000000000LL && extreme.den == 1);
    assert((extreme.cycle.size() == 4));

    return 0;
}
