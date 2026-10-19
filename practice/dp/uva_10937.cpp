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

int h, w;
vector<string> grid;

bool inside(int i, int j){
    return i >= 0 && j >= 0 && i < h && j < w;
}

/// Water, trees and every square touching a native camp, diagonals included, are off limits
bool blocked(int i, int j){
    if (grid[i][j] == '#' || grid[i][j] == '~' || grid[i][j] == '*') return true;
    for (int di = -1; di <= 1; di++){
        for (int dj = -1; dj <= 1; dj++){
            if (inside(i + di, j + dj) && grid[i + di][j + dj] == '*') return true;
        }
    }
    return false;
}

vector<vector<int>> bfs(int si, int sj, const vector<vector<bool>>& bad){
    const int dx[] = {0, 0, -1, 1}, dy[] = {1, -1, 0, 0};
    vector<vector<int>> dis(h, vector<int>(w, -1));
    queue<pair<int, int>> q;
    dis[si][sj] = 0, q.push({si, sj});
    while (!q.empty()){
        int x = q.front().first, y = q.front().second;
        q.pop();
        for (int d = 0; d < 4; d++){
            int i = x + dx[d], j = y + dy[d];
            if (inside(i, j) && !bad[i][j] && dis[i][j] == -1){
                dis[i][j] = dis[x][y] + 1;
                q.push({i, j});
            }
        }
    }
    return dis;
}

int solve(){
    vector<vector<bool>> bad(h, vector<bool>(w));
    vector<pair<int, int>> spots(1, {-1, -1});
    for (int i = 0; i < h; i++){
        for (int j = 0; j < w; j++){
            bad[i][j] = blocked(i, j);
            if (grid[i][j] != '@' && grid[i][j] != '!') continue;
            if (bad[i][j]) return -1;
            if (grid[i][j] == '@') spots[0] = {i, j};
            else spots.push_back({i, j});
        }
    }
    if (spots[0].first == -1) return -1;

    int n = spots.size();
    HamiltonianGraph g(n);
    for (int u = 0; u < n; u++){
        auto dis = bfs(spots[u].first, spots[u].second, bad);
        for (int v = 0; v < n; v++){
            int d = dis[spots[v].first][spots[v].second];
            if (d == -1) return -1;
            g.add_directed_edge(u, v, d);
        }
    }
    return g.shortest_cycle().first;
}

int main(){
    while (scanf("%d %d", &h, &w) == 2 && (h || w)){
        grid.assign(h, "");
        char row[64];
        for (int i = 0; i < h; i++){
            scanf("%s", row);
            grid[i] = row;
        }
        printf("%d\n", solve());
    }
    return 0;
}
