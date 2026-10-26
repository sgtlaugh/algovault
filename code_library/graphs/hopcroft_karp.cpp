/***
 *
 * Hopcroft-Karp
 * Maximum matching in a bipartite graph
 *
 * Complexity: O(m sqrt(n))
 *
 * HopcroftKarp hk(n_left, n_right): left vertices 0 to n_left - 1, right vertices 0 to n_right - 1
 * hk.add_edge(u, v): edge from left u to right v
 * hk.max_matching(): the size of a maximum matching, callable again after adding edges
 * hk.match_left[u]: the right vertex matched to u or -1
 * hk.match_right[v]: the left vertex matched to v or -1
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

    /// Layers left vertices by distance from free left ones, stopping at the first layer reaching a free right vertex
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
};
/// END SHARED hopcroft_karp

int main(){
    /// Left 1 only likes right 0, so left 0 must take right 1 and left 2 right 2
    HopcroftKarp hk(3, 3);
    hk.add_edge(0, 0), hk.add_edge(0, 1);
    hk.add_edge(1, 0);
    hk.add_edge(2, 1), hk.add_edge(2, 2);
    assert(hk.max_matching() == 3);
    assert(hk.match_left[1] == 0);
    assert(hk.match_right[2] == 2);

    /// More edges can be added and the matching grown
    HopcroftKarp grow(2, 2);
    grow.add_edge(0, 0), grow.add_edge(1, 0);
    assert(grow.max_matching() == 1);
    grow.add_edge(0, 1);
    assert(grow.max_matching() == 2);
    return 0;
}
