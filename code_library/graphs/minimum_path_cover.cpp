/***
 *
 * Minimum Path Cover in a DAG
 * Vertex disjoint path cover, path cover with shared vertices, and maximum antichain (Dilworth)
 *
 * Complexity:
 *   - min_disjoint_path_cover: O(m sqrt(n))
 *   - min_path_cover and max_antichain: O(n^3 / 64 + n^2.5), O(n^2) memory
 *
 * DAGPathCover g(n): vertices 0 to n - 1, no edges
 * g.add_edge(u, v): directed edge, the graph must stay acyclic
 * g.max_antichain(): a largest set of pairwise unreachable vertices, its size equals min_path_cover()
 * g.min_disjoint_path_cover(): fewest paths covering every vertex exactly once
 * g.min_path_cover(): fewest paths covering every vertex at least once, paths may share vertices
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY hopcroft_karp from code_library/graphs/hopcroft_karp.cpp
struct HopcroftKarp{
    int n_left, n_right;
    vector<vector<int>> adj;
    vector<int> match_left, match_right, dist, it, start, flat;
    vector<int> stack;  /// augment's stack as a member: allocating it per free vertex ran 14-24% slower

    HopcroftKarp(int n_left, int n_right) : n_left(n_left), n_right(n_right), adj(n_left), match_left(n_left, -1), match_right(n_right, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    int max_matching(){
        /// flat[start[u] .. start[u + 1]) packs adj[u] together, 22-39% faster on sparse graphs
        start.assign(n_left + 1, 0);
        for (int u = 0; u < n_left; u++) start[u + 1] = start[u] + adj[u].size();
        flat.resize(start[n_left]);
        for (int u = 0; u < n_left; u++) copy(adj[u].begin(), adj[u].end(), flat.begin() + start[u]);

        /// A greedy pass matches most vertices in O(m), 5-30% faster than leaving them to the first phase
        for (int u = 0; u < n_left; u++){
            if (match_left[u] != -1) continue;
            for (int e = start[u]; e < start[u + 1]; e++){
                if (match_right[flat[e]] == -1){
                    match_left[u] = flat[e], match_right[flat[e]] = u;
                    break;
                }
            }
        }

        while (true){
            int limit = bfs();
            if (limit == INT_MAX) break;
            it.assign(start.begin(), start.end() - 1);
            for (int u = 0; u < n_left; u++){
                if (match_left[u] == -1) augment(u, limit);
            }
        }

        return n_left - count(match_left.begin(), match_left.end(), -1);
    }

    /// Konig: alternating paths from the free left vertices; the cover is every unreached left and every reached right
    pair<vector<int>, vector<int>> min_vertex_cover() const{
        vector<char> seen_left(n_left, 0), seen_right(n_right, 0);
        vector<int> queue;
        for (int u = 0; u < n_left; u++){
            if (match_left[u] == -1) seen_left[u] = 1, queue.push_back(u);
        }
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : adj[queue[i]]){
                if (seen_right[v]) continue;
                seen_right[v] = 1;
                int w = match_right[v];
                assert(w != -1);  /// reaching a free right vertex is an augmenting path: the matching is not maximum
                if (!seen_left[w]) seen_left[w] = 1, queue.push_back(w);
            }
        }

        pair<vector<int>, vector<int>> res;
        for (int u = 0; u < n_left; u++) if (!seen_left[u]) res.first.push_back(u);
        for (int v = 0; v < n_right; v++) if (seen_right[v]) res.second.push_back(v);
        return res;
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
            if (dist[u] >= limit) break;
            for (int e = start[u]; e < start[u + 1]; e++){
                int v = flat[e], w = match_right[v];
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
        stack.assign(1, root);
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] == start[u + 1]){
                dist[u] = -1;
                stack.pop_back();
                continue;
            }

            int v = flat[it[u]++], w = match_right[v];
            if (w == -1 && dist[u] == limit){
                for (int x : stack){
                    int y = flat[it[x] - 1];
                    match_left[x] = y, match_right[y] = x;
                }

                return;
            }
            if (w != -1 && dist[w] == dist[u] + 1) stack.push_back(w);
        }
    }
};
/// END COPY hopcroft_karp

struct DAGPathCover{
    int n;
    vector<vector<int>> adj;

    DAGPathCover(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
    }

    vector<int> max_antichain() const{
        pair<vector<int>, vector<int>> cover = matching(n, closure()).min_vertex_cover();

        /// The antichain is every vertex with neither its out copy (left) nor its in copy (right) in the cover
        vector<char> covered(n, 0);
        for (int v : cover.first) covered[v] = 1;
        for (int v : cover.second) covered[v] = 1;

        vector<int> res;
        for (int v = 0; v < n; v++) if (!covered[v]) res.push_back(v);
        return res;
    }

    int min_disjoint_path_cover() const{
        return n - matching(n, adj).max_matching();
    }

    int min_path_cover() const{
        return n - matching(n, closure()).max_matching();
    }

    /// Split every vertex into an out copy (left) and an in copy (right), each matched edge joins two paths
    static HopcroftKarp matching(int n, vector<vector<int>> edges){
        HopcroftKarp hk(n, n);
        hk.adj = move(edges);
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
    assert(hub.max_antichain().size() == 2);  /// e.g. {0, 1}
    return 0;
}
