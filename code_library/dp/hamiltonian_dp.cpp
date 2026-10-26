/***
 *
 * Shortest Hamiltonian Path and Cycle with Bitmask DP
 * Minimum weight path or cycle visiting every vertex exactly once, with the vertex order
 *
 * Complexity: O(2^n n^2) time, O(2^n n) memory (at n = 20: 160 MiB for a path, 76 MiB for a cycle)
 *
 * HamiltonianGraph g(n): vertices 0 to n - 1, n >= 1
 * g.add_directed_edge(u, v, w): edge u -> v, parallel edges keep the minimum weight
 * g.add_edge(u, v, w): undirected edge, both directions
 * g.shortest_cycle(): {cost, order} of the cheapest cycle, order starts at 0 and closes back to it
 * g.shortest_path(): {cost, order} of the cheapest path, any start and end
 *
 * Notes:
 *   - Both return {LLONG_MAX, {}} when no such path or cycle exists
 *   - A cycle on n = 1 vertex is the self-loop 0 -> 0, on n = 2 it is 0 -> 1 -> 0
 *   - Weights may be negative and every path sum fits in a long long, so |w| <= 4e17 is safe up to n = 20
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

    pair<long long, vector<int>> solve(bool cycle) const{
        /// A cycle starts at 0, so its masks cover only vertices first..n - 1 and local vertex i is vertex i + first
        int first = cycle, k = n - first;
        if (k == 0) return w[0][0] == INF ? make_pair(INF, vector<int>{}) : make_pair(w[0][0], vector<int>{0});
        vector<long long> dp = build_dp(first);
        size_t full = (size_t(1) << k) - 1;

        long long cost = INF;
        int last = -1;
        for (int i = 0; i < k; i++){
            long long d = dp[full * k + i];
            if (d == INF || (cycle && w[i + first][0] == INF)) continue;

            long long total = cycle ? d + w[i + first][0] : d;
            if (total < cost) cost = total, last = i;
        }
        if (last == -1) return {INF, {}};

        vector<int> order(n, 0);
        size_t mask = full;
        for (int pos = n - 1; pos > first; pos--){
            order[pos] = last + first;
            size_t prev = mask ^ (size_t(1) << last);
            int from = -1;
            for (int j = 0; j < k && from == -1; j++){
                long long d = dp[prev * k + j], step = w[j + first][last + first];
                if (d != INF && step != INF && d + step == dp[mask * k + last]) from = j;
            }
            mask = prev, last = from;
        }
        order[first] = last + first;

        return {cost, order};
    }

    vector<long long> build_dp(int first) const{
        int k = n - first;
        vector<long long> dp((size_t(1) << k) * k, INF);
        vector<size_t> into(k, 0);
        vector<vector<long long>> w_into(k, vector<long long>(k));
        for (int j = 0; j < k; j++){
            for (int i = 0; i < k; i++){
                w_into[i][j] = w[j + first][i + first];
                if (w_into[i][j] != INF) into[i] |= size_t(1) << j;
            }
        }

        for (int i = 0; i < k; i++) dp[(size_t(1) << i) * k + i] = first ? w[0][i + first] : 0;

        for (size_t mask = 1; mask < (size_t(1) << k); mask++){
            for (int i = 0; i < k; i++){
                if (!(mask >> i & 1) || mask == (size_t(1) << i)) continue;

                size_t prev = mask ^ (size_t(1) << i);
                const vector<long long>& row = w_into[i];
                long long best = INF;
                for (size_t rest = prev & into[i]; rest; rest &= rest - 1){
                    int j = __builtin_ctzll(rest);
                    long long d = dp[prev * k + j];
                    if (d != INF) best = min(best, d + row[j]);
                }
                dp[mask * k + i] = best;
            }
        }

        return dp;
    }
};

int main(){
    HamiltonianGraph g(3);
    vector<vector<long long>> mat = {{8, 1, 6}, {3, 5, 7}, {4, 9, 2}};
    for (int u = 0; u < 3; u++){
        for (int v = 0; v < 3; v++) g.add_directed_edge(u, v, mat[u][v]);
    }
    assert((g.shortest_path() == pair<long long, vector<int>>{5, {2, 0, 1}}));    /// 4 + 1
    assert((g.shortest_cycle() == pair<long long, vector<int>>{12, {0, 1, 2}}));  /// 1 + 7 + 4

    HamiltonianGraph chain(3);
    chain.add_directed_edge(0, 1, 2), chain.add_directed_edge(1, 2, 3);
    assert(chain.shortest_cycle().first == LLONG_MAX);  /// no edge back to 0
    chain.add_directed_edge(2, 0, -9);
    assert((chain.shortest_cycle() == pair<long long, vector<int>>{-4, {0, 1, 2}}));  /// negative weights are fine
    return 0;
}
