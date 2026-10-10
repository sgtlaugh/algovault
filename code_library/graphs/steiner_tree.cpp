/***
 *
 * Steiner Tree (Dreyfus Wagner)
 * Minimum total weight of edges connecting every terminal of an undirected graph
 *
 * Complexity: O(3^k n + 2^k m log n) time, O(2^k n) memory for k terminals
 *
 * SteinerTree g(n): 0-based vertices
 * g.add_edge(u, v, w): undirected edge with w >= 0
 * g.solve(terminals): the minimum weight, 0 for fewer than two terminals, STEINER_INF if they are not all connected
 *     does not modify the graph, so one graph answers several terminal sets
 *     the total weight of all edges must stay below STEINER_INF (~2.3e18)
 *
 * dp[mask][v] = cheapest tree connecting the terminals in mask plus vertex v:
 * merge two smaller trees at v, then spread outward with a multi-source Dijkstra
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long STEINER_INF = numeric_limits<long long>::max() / 4;

struct SteinerTree{
    int n;
    vector<vector<pair<int, long long>>> adj;

    SteinerTree(int n) : n(n), adj(n) {}

    void add_edge(int u, int v, long long w){
        assert(w >= 0);
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    long long solve(const vector<int>& terminals) const{
        int k = terminals.size();
        if (k <= 1) return 0;

        vector<vector<long long>> dp(1 << k, vector<long long>(n, STEINER_INF));
        for (int i = 0; i < k; i++) dp[1 << i][terminals[i]] = 0;

        for (int mask = 1; mask < (1 << k); mask++){
            vector<long long>& d = dp[mask];
            for (int sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask){
                if (sub < (sub ^ mask)) continue;  /// each split once
                for (int v = 0; v < n; v++) d[v] = min(d[v], dp[sub][v] + dp[sub ^ mask][v]);
            }

            priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
            for (int v = 0; v < n; v++){
                if (d[v] < STEINER_INF) heap.push({d[v], v});
            }
            while (!heap.empty()){
                auto [dv, v] = heap.top();
                heap.pop();
                if (dv != d[v]) continue;
                for (auto [u, w] : adj[v]){
                    if (dv + w < d[u]) d[u] = dv + w, heap.push({d[u], u});
                }
            }
        }

        return *min_element(dp[(1 << k) - 1].begin(), dp[(1 << k) - 1].end());
    }
};

int main(){
    /***
     * Star with a cheap hub: terminals 0, 1, 2 joined through vertex 3 at cost 1 each
     * The triangle edges among terminals cost 5, so the best tree uses the hub, total 3
    ***/
    SteinerTree star(4);
    star.add_edge(0, 3, 1), star.add_edge(1, 3, 1), star.add_edge(2, 3, 1);
    star.add_edge(0, 1, 5), star.add_edge(1, 2, 5), star.add_edge(0, 2, 5);
    assert(star.solve({0, 1, 2}) == 3);
    assert(star.solve({0, 1}) == 2);  /// the same graph answers another terminal set

    SteinerTree path(5);
    path.add_edge(0, 1, 4), path.add_edge(1, 2, 6), path.add_edge(2, 3, 1);
    assert(path.solve({0, 3}) == 11);
    assert(path.solve({0, 4}) == STEINER_INF);  /// 4 is isolated
    return 0;
}
