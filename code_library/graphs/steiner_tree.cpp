/***
 *
 * Steiner Tree (Dreyfus Wagner)
 * Minimum total weight of edges connecting every terminal of an undirected graph
 *
 * Complexity: O(3^k n + 2^k m log n) time, O(2^k n) memory for k terminals
 *
 * steiner_tree(n, edges, terminals): edges[i] = {u, v, w} with w >= 0, 0-based vertices
 *     returns the minimum weight, 0 for fewer than two terminals, STEINER_INF if they are not all connected
 *     the total weight of all edges must stay below STEINER_INF (~2.3e18)
 *
 * dp[mask][v] = cheapest tree connecting the terminals in mask plus vertex v:
 * merge two smaller trees at v, then spread outward with a multi-source Dijkstra
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long STEINER_INF = numeric_limits<long long>::max() / 4;

long long steiner_tree(int n, const vector<array<long long, 3>>& edges, const vector<int>& terminals){
    int k = terminals.size();
    if (k <= 1) return 0;

    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges){
        assert(w >= 0);
        adj[u].push_back({(int)v, w});
        adj[v].push_back({(int)u, w});
    }

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

int main(){
    /***
     * Star with a cheap hub: terminals 0, 1, 2 joined through vertex 3 at cost 1 each
     * The triangle edges among terminals cost 5, so the best tree uses the hub, total 3
    ***/
    vector<array<long long, 3>> star = {{0, 3, 1}, {1, 3, 1}, {2, 3, 1}, {0, 1, 5}, {1, 2, 5}, {0, 2, 5}};
    assert(steiner_tree(4, star, {0, 1, 2}) == 3);
    assert(steiner_tree(4, star, {0, 1}) == 2);
    assert(steiner_tree(4, star, {2}) == 0);
    assert(steiner_tree(4, star, {}) == 0);

    vector<array<long long, 3>> path = {{0, 1, 4}, {1, 2, 6}, {2, 3, 1}};
    assert(steiner_tree(4, path, {0, 3}) == 11);
    assert(steiner_tree(5, path, {0, 4}) == STEINER_INF);
    assert(steiner_tree(4, path, {3, 2, 3}) == 1);
    return 0;
}
