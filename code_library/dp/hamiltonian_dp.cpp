/***
 *
 * Shortest Hamiltonian Path and Cycle with Bitmask DP
 * Minimum weight path or cycle visiting every vertex exactly once, with the vertex order
 *
 * Complexity: O(2^n * n^2) time, O(2^n * n) memory (2^n * n long longs, 160 MiB at n = 20)
 *
 * HamiltonianGraph g(n): vertices 0 to n - 1, n >= 1
 * g.add_directed_edge(u, v, w): edge u -> v, parallel edges keep the minimum weight
 * g.add_edge(u, v, w): undirected edge, both directions
 * g.shortest_path(): {cost, order} of the cheapest path, any start and end
 * g.shortest_cycle(): {cost, order} of the cheapest cycle, order starts at 0 and closes back to it
 *     a cycle on n = 1 vertex is the self-loop 0 -> 0, on n = 2 it is 0 -> 1 -> 0
 * Both return {LLONG_MAX, {}} when no such path or cycle exists
 *
 * Weights may be negative; every path sum must fit in a long long, so |w| <= 4e17 is safe up to n = 20
 * Missing edges are skipped instead of stored as a large weight, so no INF + w overflow
 *
 * dp[mask][i] = cheapest path visiting exactly the vertices in mask and ending at i
 * The cycle runs the same DP with only dp[{0}][0] seeded, then closes with the edge back to 0
 *
 * Example:
 *   HamiltonianGraph g(3);
 *   g.add_directed_edge(2, 0, 4), g.add_directed_edge(0, 1, 1);
 *   auto [cost, order] = g.shortest_path();  // cost = 5, order = {2, 0, 1}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct HamiltonianGraph{
    static constexpr long long INF = LLONG_MAX;

    int n;
    vector<vector<long long>> w;

    HamiltonianGraph(int n) : n(n), w(n, vector<long long>(n, INF)) {}

    void add_directed_edge(int u, int v, long long weight){
        w[u][v] = min(w[u][v], weight);
    }

    void add_edge(int u, int v, long long weight){
        add_directed_edge(u, v, weight);
        add_directed_edge(v, u, weight);
    }

    pair<long long, vector<int>> shortest_cycle() const{
        return solve(true);
    }

    pair<long long, vector<int>> shortest_path() const{
        return solve(false);
    }

private:
    vector<long long> build_dp(bool cycle) const{
        vector<long long> dp((size_t(1) << n) * n, INF);
        if (cycle) dp[n] = 0;  // flat index mask * n + i with mask = {0}, i = 0
        else{
            for (int i = 0; i < n; i++) dp[(size_t(1) << i) * n + i] = 0;
        }

        for (size_t mask = 1; mask < (size_t(1) << n); mask++){
            for (int i = 0; i < n; i++){
                if (!(mask >> i & 1) || mask == (size_t(1) << i)) continue;

                size_t prev = mask ^ (size_t(1) << i);
                long long best = INF;
                for (int j = 0; j < n; j++){
                    long long d = dp[prev * n + j];
                    if (d != INF && w[j][i] != INF) best = min(best, d + w[j][i]);
                }
                dp[mask * n + i] = best;
            }
        }

        return dp;
    }

    pair<long long, vector<int>> solve(bool cycle) const{
        vector<long long> dp = build_dp(cycle);
        size_t full = (size_t(1) << n) - 1;

        long long cost = INF;
        int last = -1;
        for (int i = 0; i < n; i++){
            long long d = dp[full * n + i];
            if (d == INF || (cycle && w[i][0] == INF)) continue;

            long long total = cycle ? d + w[i][0] : d;
            if (total < cost) cost = total, last = i;
        }
        if (last == -1) return {INF, {}};

        vector<int> order(n);
        size_t mask = full;
        for (int pos = n - 1; pos > 0; pos--){
            order[pos] = last;
            size_t prev = mask ^ (size_t(1) << last);
            int from = -1;
            for (int j = 0; j < n && from == -1; j++){
                long long d = dp[prev * n + j];
                if (d != INF && w[j][last] != INF && d + w[j][last] == dp[mask * n + last]) from = j;
            }
            mask = prev, last = from;
        }
        order[0] = last;

        return {cost, order};
    }
};

int main(){
    HamiltonianGraph a(3);
    vector<vector<long long>> mat = {{8, 1, 6}, {3, 5, 7}, {4, 9, 2}};
    for (int u = 0; u < 3; u++){
        for (int v = 0; v < 3; v++) a.add_directed_edge(u, v, mat[u][v]);
    }
    assert((a.shortest_path() == pair<long long, vector<int>>{5, {2, 0, 1}}));
    assert((a.shortest_cycle() == pair<long long, vector<int>>{12, {0, 1, 2}}));

    HamiltonianGraph b(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 3}, {1, 4}, {2, 3}, {2, 4}}) b.add_edge(u, v, 1);
    for (int u = 0; u < 5; u++){
        for (int v = u + 1; v < 5; v++) b.add_edge(u, v, 10);
    }
    auto [cost, order] = b.shortest_cycle();
    assert(cost == 5);
    assert((order == vector<int>{0, 1, 4, 2, 3} || order == vector<int>{0, 3, 2, 4, 1}));

    HamiltonianGraph single(1);
    assert((single.shortest_path() == pair<long long, vector<int>>{0, {0}}));
    assert(single.shortest_cycle().second.empty());
    single.add_directed_edge(0, 0, 7);
    assert((single.shortest_cycle() == pair<long long, vector<int>>{7, {0}}));

    HamiltonianGraph chain(3);
    chain.add_directed_edge(0, 1, 2), chain.add_directed_edge(1, 2, 3);
    assert((chain.shortest_path() == pair<long long, vector<int>>{5, {0, 1, 2}}));
    assert((chain.shortest_cycle() == pair<long long, vector<int>>{LLONG_MAX, {}}));
    chain.add_directed_edge(2, 0, -9);
    assert((chain.shortest_cycle() == pair<long long, vector<int>>{-4, {0, 1, 2}}));
    assert((chain.shortest_path() == pair<long long, vector<int>>{-7, {2, 0, 1}}));

    HamiltonianGraph split(4);
    split.add_edge(0, 1, 1), split.add_edge(2, 3, 1);
    assert((split.shortest_path() == pair<long long, vector<int>>{LLONG_MAX, {}}));

    HamiltonianGraph pair_graph(2);
    pair_graph.add_edge(0, 1, 6), pair_graph.add_directed_edge(0, 1, 4);
    assert((pair_graph.shortest_path() == pair<long long, vector<int>>{4, {0, 1}}));
    assert((pair_graph.shortest_cycle() == pair<long long, vector<int>>{10, {0, 1}}));

    return 0;
}
