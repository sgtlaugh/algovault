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

constexpr long long HamiltonianGraph::INF;

const int MAXD = 1000, OFF = 4, SZ = MAXD + 2 * OFF;
int dist_inf[SZ][SZ];

void knight_bfs(){
    const int dx[] = {1, 1, -1, -1, 2, 2, -2, -2}, dy[] = {2, -2, 2, -2, 1, -1, 1, -1};

    memset(dist_inf, -1, sizeof(dist_inf));
    queue<pair<int, int>> q;
    dist_inf[OFF][OFF] = 0, q.push(make_pair(OFF, OFF));
    while (!q.empty()){
        int x = q.front().first, y = q.front().second;
        q.pop();
        for (int k = 0; k < 8; k++){
            int a = x + dx[k], b = y + dy[k];
            if (a >= 0 && b >= 0 && a < SZ && b < SZ && dist_inf[a][b] == -1){
                dist_inf[a][b] = dist_inf[x][y] + 1, q.push(make_pair(a, b));
            }
        }
    }
}

int knight_dist(int n, int r1, int c1, int r2, int c2){
    auto corner = [n](int r, int c){ return (r == 1 || r == n) && (c == 1 || c == n); };
    int dr = abs(r1 - r2), dc = abs(c1 - c2);
    if (dr == 1 && dc == 1 && (corner(r1, c1) || corner(r2, c2))) return 4;
    if (n == 4 && corner(r1, c1) && corner(r2, c2) && dr + dc == 3) return 5;
    return dist_inf[OFF + dr][OFF + dc];
}

int main(){
    knight_bfs();

    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n, k;
        scanf("%d %d", &n, &k);
        vector<int> r(k), c(k);
        for (int i = 0; i < k; i++) scanf("%d %d", &r[i], &c[i]);

        HamiltonianGraph g(k);
        for (int i = 0; i < k; i++){
            for (int j = 0; j < k; j++) g.add_directed_edge(i, j, knight_dist(n, r[i], c[i], r[j], c[j]));
        }
        printf("Case %d: %lld\n", tc, g.shortest_cycle().first);
    }
    return 0;
}
