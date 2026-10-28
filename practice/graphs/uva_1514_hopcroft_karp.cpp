/***
 *
 * Hopcroft Karp
 * Maximum matching in a bipartite graph
 *
 * Complexity: O(m * sqrt(n)), O(n + m) memory
 *
 * HopcroftKarp hk(n_left, n_right): left vertices 0..n_left - 1, right vertices 0..n_right - 1
 * hk.add_edge(u, v): edge from left u to right v
 * hk.max_matching(): size of a maximum matching, can be called again after adding more edges
 * hk.match_left[u]: right vertex matched to u or -1, hk.match_right[v]: left vertex matched to v or -1
 *
 * Each phase augments along shortest paths only, which is what bounds the phases by O(sqrt(n))
 * The augmenting search is iterative, so long augmenting paths cannot overflow the stack
 *
 * minimum_path_cover.cpp carries a checked copy of this struct
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED hopcroft_karp
struct HopcroftKarp{
    int n_left, n_right;
    vector<vector<int>> adj;
    vector<int> match_left, match_right, dist, it;

    HopcroftKarp(int n_left, int n_right) : n_left(n_left), n_right(n_right), adj(n_left), match_left(n_left, -1), match_right(n_right, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    /// Layers left vertices by distance from free left vertices, stopping at the first layer that reaches a free right vertex
    int bfs(){
        dist.assign(n_left, -1);
        vector<int> queue;
        for (int u = 0; u < n_left; u++){
            if (match_left[u] == -1) dist[u] = 0, queue.push_back(u);
        }

        int limit = INT_MAX;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int v : adj[u]){
                int w = match_right[v];
                if (w == -1) limit = min(limit, dist[u]);
                else if (dist[w] == -1 && dist[u] < limit){
                    dist[w] = dist[u] + 1;
                    queue.push_back(w);
                }
            }
        }

        return limit;
    }

    void augment(int root, int limit){
        vector<int> stack = {root};
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == (int)adj[u].size()){
                dist[u] = -1;
                stack.pop_back();
                continue;
            }

            int v = adj[u][it[u]++], w = match_right[v];
            if (w == -1 && dist[u] == limit){
                for (int x : stack){
                    int y = adj[x][it[x] - 1];
                    match_left[x] = y, match_right[y] = x;
                }

                return;
            }
            if (w != -1 && dist[w] == dist[u] + 1) stack.push_back(w);
        }
    }

    int max_matching(){
        while (true){
            int limit = bfs();
            if (limit == INT_MAX) break;
            it.assign(n_left, 0);
            for (int u = 0; u < n_left; u++){
                if (match_left[u] == -1) augment(u, limit);
            }
        }

        return n_left - count(match_left.begin(), match_left.end(), -1);
    }
};
/// END SHARED hopcroft_karp

int main(){
    int t;

    scanf("%d", &t);
    while (t--){
        int n, m;
        scanf("%d %d", &n, &m);
        vector<string> grid(n);
        char row[512];
        for (auto& s: grid) scanf("%s", row), s = row;

        int blacks = 0, whites = 0;
        vector<vector<int>> id(n, vector<int>(m, -1));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                if (grid[i][j] == 'B') id[i][j] = blacks++;
                if (grid[i][j] == 'W') id[i][j] = whites++;
            }
        }
        if (2 * blacks != whites){
            puts("NO");
            continue;
        }

        HopcroftKarp hk(2 * blacks, whites);
        const int dx[] = {0, 0, -1, 1}, dy[] = {-1, 1, 0, 0};
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                if (grid[i][j] != 'B') continue;
                for (int k = 0; k < 4; k++){
                    int x = i + dx[k], y = j + dy[k];
                    if (x >= 0 && x < n && y >= 0 && y < m && grid[x][y] == 'W'){
                        hk.add_edge(id[i][j] + (k < 2 ? 0 : blacks), id[x][y]);
                    }
                }
            }
        }
        puts(hk.max_matching() == whites ? "YES" : "NO");
    }

    return 0;
}
