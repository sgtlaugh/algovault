/***
 *
 * Shortest Hamiltonian Path and Cycle with Bitmask DP
 * Minimum weight path or cycle visiting every vertex exactly once, with the vertex order
 *
 * Complexity: O(2^n n^2) time, O(2^n n) memory (160 MiB at n = 20)
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

    vector<long long> build_dp(bool cycle) const{
        vector<long long> dp((size_t(1) << n) * n, INF);
        vector<size_t> into(n, 0);
        vector<vector<long long>> w_into(n, vector<long long>(n));
        for (int j = 0; j < n; j++){
            for (int i = 0; i < n; i++){
                w_into[i][j] = w[j][i];
                if (w[j][i] != INF) into[i] |= size_t(1) << j;
            }
        }

        if (cycle) dp[n] = 0;  // flat index mask * n + i with mask = {0}, i = 0
        else{
            for (int i = 0; i < n; i++) dp[(size_t(1) << i) * n + i] = 0;
        }

        /// A cycle starts at vertex 0, so every finite cycle state contains it and the other masks are skipped
        for (size_t mask = 1; mask < (size_t(1) << n); mask++){
            if (cycle && !(mask & 1)) continue;
            for (int i = 0; i < n; i++){
                if (!(mask >> i & 1) || mask == (size_t(1) << i)) continue;

                size_t prev = mask ^ (size_t(1) << i);
                const vector<long long>& row = w_into[i];
                long long best = INF;
                for (size_t rest = prev & into[i]; rest; rest &= rest - 1){
                    int j = __builtin_ctzll(rest);
                    long long d = dp[prev * n + j];
                    if (d != INF) best = min(best, d + row[j]);
                }
                dp[mask * n + i] = best;
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
