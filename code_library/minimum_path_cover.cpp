/***
 *
 * Minimum Path Cover in a DAG
 * Vertex disjoint path cover, path cover with shared vertices, and maximum antichain (Dilworth)
 *
 * Complexity:
 *   - min_disjoint_path_cover: O(m sqrt(n)) via bipartite matching
 *   - min_path_cover and max_antichain: O(n^3 / 64) transitive closure, then matching on up to n^2 edges
 *
 * DAGPathCover g(n); g.add_edge(u, v): directed edge, the graph must be acyclic
 * g.min_disjoint_path_cover(): fewest paths covering every vertex exactly once
 * g.min_path_cover(): fewest paths covering every vertex at least once, paths may share vertices
 * g.max_antichain(): largest set of pairwise unreachable vertices, its size equals min_path_cover() by Dilworth's theorem
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Same as hopcroft_karp.cpp, shortest augmenting paths per phase and an iterative search
struct HopcroftKarp{
    int n_left, n_right;
    vector<vector<int>> adj;
    vector<int> match_left, match_right, dist, it;

    HopcroftKarp(int n_left, int n_right) : n_left(n_left), n_right(n_right), adj(n_left), match_left(n_left, -1), match_right(n_right, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

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
                else if (dist[w] == -1 && dist[u] < limit) dist[w] = dist[u] + 1, queue.push_back(w);
            }
        }
        return limit;
    }

    bool augment(int root, int limit){
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
                return true;
            }
            if (w != -1 && dist[w] == dist[u] + 1) stack.push_back(w);
        }
        return false;
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

struct DAGPathCover{
    int n;
    vector<vector<int>> adj;

    DAGPathCover(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    /// Split every vertex into an out copy (left) and an in copy (right), each matched edge joins two paths
    static HopcroftKarp matching(int n, const vector<vector<int>>& edges){
        HopcroftKarp hk(n, n);
        for (int u = 0; u < n; u++){
            for (int v : edges[u]) hk.add_edge(u, v);
        }
        hk.max_matching();
        return hk;
    }

    vector<vector<int>> closure() const{
        int words = (n + 63) / 64;
        vector<vector<unsigned long long>> reach(n, vector<unsigned long long>(words, 0));

        /// Reverse topological order, so every successor's reach set is complete before it is merged
        vector<int> indegree(n, 0), order;
        for (int u = 0; u < n; u++){
            for (int v : adj[u]) indegree[v]++;
        }
        for (int u = 0; u < n; u++){
            if (!indegree[u]) order.push_back(u);
        }
        for (int i = 0; i < (int)order.size(); i++){
            for (int v : adj[order[i]]){
                if (!--indegree[v]) order.push_back(v);
            }
        }
        assert((int)order.size() == n);

        for (int i = n - 1; i >= 0; i--){
            int u = order[i];
            for (int v : adj[u]){
                reach[u][v >> 6] |= 1ULL << (v & 63);
                for (int w = 0; w < words; w++) reach[u][w] |= reach[v][w];
            }
        }

        vector<vector<int>> res(n);
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++){
                if (reach[u][v >> 6] >> (v & 63) & 1) res[u].push_back(v);
            }
        }
        return res;
    }

    int min_disjoint_path_cover() const{
        return n - matching(n, adj).max_matching();
    }

    int min_path_cover() const{
        return n - matching(n, closure()).max_matching();
    }

    vector<int> max_antichain() const{
        vector<vector<int>> edges = closure();
        HopcroftKarp hk = matching(n, edges);

        /// Konig: alternating reachability from free left vertices gives the minimum vertex cover
        vector<char> seen_left(n, 0), seen_right(n, 0);
        vector<int> queue;
        for (int u = 0; u < n; u++){
            if (hk.match_left[u] == -1) seen_left[u] = 1, queue.push_back(u);
        }
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : edges[queue[i]]){
                if (seen_right[v]) continue;
                seen_right[v] = 1;
                int w = hk.match_right[v];
                if (w != -1 && !seen_left[w]) seen_left[w] = 1, queue.push_back(w);
            }
        }

        /// The cover is unvisited left plus visited right copies, the antichain is every vertex with neither copy in it
        vector<int> res;
        for (int v = 0; v < n; v++){
            if (seen_left[v] && !seen_right[v]) res.push_back(v);
        }
        return res;
    }
};

int main(){
    /***
     * 0 -> 1 -> 3 and 0 -> 2 -> 3, plus a separate vertex 4
     * Disjoint cover: 0-1-3, 2, 4 = 3 paths; shared cover: 0-1-3, 0-2-3, 4 = 3 paths; antichain {1, 2, 4}
    ***/
    DAGPathCover g(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 3}, {0, 2}, {2, 3}}) g.add_edge(u, v);
    assert(g.min_disjoint_path_cover() == 3);
    assert(g.min_path_cover() == 3);
    assert((g.max_antichain() == vector<int>{1, 2, 4}));

    /***
     * Star through a hub: 0 -> 2, 1 -> 2, 2 -> 3, 2 -> 4
     * Disjoint: 0-2-3, 1, 4 = 3 paths; shared vertices allowed: 0-2-3, 1-2-4 = 2 paths
    ***/
    DAGPathCover hub(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 2}, {1, 2}, {2, 3}, {2, 4}}) hub.add_edge(u, v);
    assert(hub.min_disjoint_path_cover() == 3);
    assert(hub.min_path_cover() == 2);
    assert(hub.max_antichain().size() == 2);

    DAGPathCover lonely(3);
    assert(lonely.min_disjoint_path_cover() == 3 && lonely.min_path_cover() == 3 && lonely.max_antichain().size() == 3);

    return 0;
}
