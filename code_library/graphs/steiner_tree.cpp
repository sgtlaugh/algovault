/***
 *
 * Steiner Tree (Dreyfus Wagner)
 * Cheapest tree connecting every terminal of an undirected graph, paying for its edges and its vertices
 *
 * Complexity: O(3^k n + 2^k m log n) time, O(2^k n) memory for k terminals
 *
 * SteinerTree g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v, w): adds an undirected edge with w >= 0 and returns its id, ids count up from 0 in add_edge order
 * g.solve(terminals): minimum edge weight plus vertex cost of a tree joining them, STEINER_INF if not all connected
 * g.solve_subsets(terminals): best, where best[mask] is solve of the terminals[i] with bit i set in mask
 * g.tree(terminals): distinct edge ids of a cheapest tree, the terminals must be connected
 * g.cost: cost[v] >= 0 is paid once when v is in the tree, all 0 until set
 *
 * Notes:
 *   - solve and tree do not modify the graph, so one graph answers several terminal sets
 *   - The tree's vertices are the endpoints of tree's edges, or the lone terminal when there is no edge
 *   - solve of no terminals is 0
 *   - The edges of tree form a tree when all weights are positive, zero weight edges can close a cycle
 *   - Steiner forest, each group connected: forest[s] = min over t of forest[s ^ t] + best[t], t a union of groups in s
 *   - The total weight of all edges plus the total cost must stay below STEINER_INF (~2.3e18)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long STEINER_INF = numeric_limits<long long>::max() / 4;

struct SteinerTree{
    struct Edge{
        int to, id;
        long long w;
    };

    int n, m = 0;
    vector<vector<Edge>> adj;
    vector<long long> cost;

    SteinerTree(int n) : n(n), adj(n), cost(n, 0) {}

    int add_edge(int u, int v, long long w){
        assert(w >= 0);
        adj[u].push_back({v, m, w});
        adj[v].push_back({u, m, w});
        return m++;
    }

    long long solve(const vector<int>& terminals) const{
        return solve_subsets(terminals).back();
    }

    vector<long long> solve_subsets(const vector<int>& terminals) const{
        int k = terminals.size();
        if (k == 0) return {0};

        vector<vector<long long>> dp = table(terminals);
        vector<long long> best(1 << k, STEINER_INF);
        best[0] = 0;
        for (int mask = 1; mask < (1 << k); mask++){
            for (int v = 0; v < n; v++) best[mask] = min(best[mask], dp[mask][v] + cost[v]);
        }

        return best;
    }

    vector<int> tree(const vector<int>& terminals) const{
        int k = terminals.size();
        if (k <= 1) return {};

        vector<vector<long long>> dp = table(terminals);
        int full = (1 << k) - 1;
        int root = 0;
        for (int v = 0; v < n; v++){
            if (dp[full][v] + cost[v] < dp[full][root] + cost[root]) root = v;
        }
        assert(dp[full][root] < STEINER_INF);

        vector<int> edges;
        vector<bool> taken(m, false);
        vector<pair<int, int>> todo = {{full, root}};
        while (!todo.empty()){
            auto [mask, v] = todo.back();
            todo.pop_back();
            const vector<long long>& d = dp[mask];

            /// breadth first over edges with d[u] = d[to] + cost[to] + w, a greedy walk can cycle on zero weight edges
            vector<int> parent(n, -1), parent_edge(n, -1), queue = {v};
            parent[v] = v;
            int x = -1, sub = 0;
            for (int i = 0; ; i++){
                int u = queue[i];
                for (sub = (mask - 1) & mask; sub > 0; sub = (sub - 1) & mask){
                    if (dp[sub][u] + dp[sub ^ mask][u] == d[u]) break;
                }
                if (sub > 0 || (mask == (mask & -mask) && u == terminals[__lg(mask)])){
                    x = u;
                    break;
                }

                for (auto [to, id, w] : adj[u]){
                    if (parent[to] != -1 || d[to] + cost[to] + w != d[u]) continue;
                    parent[to] = u, parent_edge[to] = id;
                    queue.push_back(to);
                }
            }

            for (int u = x; u != v; u = parent[u]){
                if (!taken[parent_edge[u]]) taken[parent_edge[u]] = true, edges.push_back(parent_edge[u]);
            }
            if (sub > 0) todo.push_back({sub, x}), todo.push_back({sub ^ mask, x});
        }

        return edges;
    }

    vector<vector<long long>> table(const vector<int>& terminals) const{
        int k = terminals.size();
        for (long long c : cost) assert(c >= 0);

        /// dp[mask][v] leaves out cost[v], so merging two trees at v pays it once without a subtraction
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

                dv += cost[v];
                for (auto [u, id, w] : adj[v]){
                    if (dv + w < d[u]) d[u] = dv + w, heap.push({d[u], u});
                }
            }
        }

        return dp;
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

    vector<int> hub = star.tree({0, 1, 2});
    sort(hub.begin(), hub.end());
    assert(hub == vector<int>({0, 1, 2}));

    SteinerTree path(5);
    path.add_edge(0, 1, 4), path.add_edge(1, 2, 6), path.add_edge(2, 3, 1);
    assert(path.solve({0, 3}) == 11);
    assert(path.solve({0, 4}) == STEINER_INF);  /// 4 is isolated

    /// Bit i of the mask picks terminals[i], so mask 6 is {1, 3} at cost 6 + 1
    vector<long long> best = path.solve_subsets({0, 1, 3});
    assert(best == vector<long long>({0, 0, 0, 4, 0, 11, 7, 11}));

    /// Free edges and paid vertices: 0 reaches 2 through 1 at 1 + 5 + 1 or through 3 at 1 + 2 + 1
    SteinerTree cells(4);
    cells.cost = {1, 5, 1, 2};
    cells.add_edge(0, 1, 0), cells.add_edge(1, 2, 0), cells.add_edge(0, 3, 0), cells.add_edge(3, 2, 0);
    assert(cells.solve({0, 2}) == 4);
    assert(cells.solve({1}) == 5);  /// a lone terminal pays its own cost

    vector<int> detour = cells.tree({0, 2});
    sort(detour.begin(), detour.end());
    assert(detour == vector<int>({2, 3}));
    return 0;
}
