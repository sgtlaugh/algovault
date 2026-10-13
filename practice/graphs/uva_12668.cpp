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
    int n;
    while (scanf("%d", &n) == 1){
        vector<string> board(n);
        char line[128];
        for (auto& row: board){
            scanf("%127s", line);
            row = line;
        }

        /// Pawns cut rows and columns into segments, a rook takes one row segment and one column segment
        vector<vector<int>> row_seg(n, vector<int>(n, -1)), col_seg(n, vector<int>(n, -1));
        int rows = 0, cols = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (board[i][j] == 'X') continue;
                row_seg[i][j] = (j && board[i][j - 1] == '.') ? row_seg[i][j - 1] : rows++;
                col_seg[i][j] = (i && board[i - 1][j] == '.') ? col_seg[i - 1][j] : cols++;
            }
        }

        HopcroftKarp hk(rows, cols);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (board[i][j] == '.') hk.add_edge(row_seg[i][j], col_seg[i][j]);
            }
        }
        printf("%d\n", hk.max_matching());
    }
    return 0;
}
