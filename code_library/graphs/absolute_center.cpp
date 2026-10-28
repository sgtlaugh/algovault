/***
 *
 * Absolute Center (Kariv Hakimi)
 * The point of a weighted undirected graph, on a vertex or inside an edge, whose farthest vertex is nearest
 *
 * Complexity: O(n m log n + n^2 log n) time, O(n^2 + m) memory
 *
 * AbsoluteCenter g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v, w): adds an undirected edge with w >= 0 and returns its id, ids count up from 0 in add_edge order
 * g.solve(): {diameter, edge, offset, tree} for the center, the graph must be connected
 *
 * Notes:
 *   - The center lies on edge id edge at distance offset / 2 from that edge's u, so offset is in [0, 2 w]
 *   - diameter is twice the radius, the center's distance to its farthest vertex, and is the diameter of tree
 *   - tree holds the n - 1 edge ids of a minimum diameter spanning tree, keeping each vertex's distance to the center
 *   - A graph without edges, so n = 1, gives {0, -1, 0, {}}
 *   - Parallel edges, self loops and zero weights are allowed
 *   - (n - 1) * max w must stay below 2^61
 *   - Method: Kariv and Hakimi, An Algorithmic Approach to Network Location Problems. I: The p-Centers, 1979
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct AbsoluteCenterResult{
    long long diameter;
    int edge;
    long long offset;
    vector<int> tree;
};

struct AbsoluteCenter{
    /// above every doubled distance, so the shortest path tree's Dijkstra needs it this large
    static constexpr long long INF = numeric_limits<long long>::max() / 2;

    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    vector<Edge> edges;
    vector<vector<pair<int, int>>> adj;

    AbsoluteCenter(int n) : n(n), adj(n) {}

    int add_edge(int u, int v, long long w){
        assert(w >= 0);
        int id = edges.size();
        edges.push_back({u, v, w});
        adj[u].push_back({v, id});
        adj[v].push_back({u, id});
        return id;
    }

    AbsoluteCenterResult solve() const{
        vector<vector<long long>> dist(n);
        for (int s = 0; s < n; s++){
            vector<long long> from(n, INF);
            from[s] = 0;
            dist[s] = dijkstra(from, 1).first;
        }
        for (int v = 0; v < n; v++) assert(dist[0][v] < INF);

        /// A point at distance t from u reaches x in min(du[x] + t, dv[x] + w - t), and the farthest x decides
        AbsoluteCenterResult res{0, -1, 0, {}};
        long long best = INF;
        vector<long long> ecc(n);
        for (int u = 0; u < n; u++) ecc[u] = *max_element(dist[u].begin(), dist[u].end());

        vector<int> order(n);
        for (int u = 0; u < n; u++){
            const vector<long long>& du = dist[u];
            bool sorted = false;

            for (auto [v, id] : adj[u]){
                if (edges[id].u != u) continue;
                const vector<long long>& dv = dist[v];
                long long w = edges[id].w;

                /// Every point of the edge has radius >= max(ecc[u] - t, ecc[v] - w + t), so twice it is >= this sum
                if (ecc[u] + ecc[v] - w >= best) continue;
                if (!sorted){
                    iota(order.begin(), order.end(), 0);
                    sort(order.begin(), order.end(), [&](int a, int b){ return du[a] > du[b]; });
                    sorted = true;
                }

                /// p has the largest dv so far, the optimum on the edge is t = 0, t = w or where p's tent meets x's
                int p = order[0];
                if (2 * du[p] < best) best = 2 * du[p], res.edge = id, res.offset = 0;
                for (int x : order){
                    if (dv[x] <= dv[p]) continue;
                    long long meet = du[x] + dv[p] + w;
                    if (meet < best) best = meet, res.edge = id, res.offset = dv[p] + w - du[x];
                    p = x;
                }
                if (2 * dv[p] < best) best = 2 * dv[p], res.edge = id, res.offset = 2 * w;
            }
        }
        if (res.edge == -1) return res;

        res.diameter = best;
        res.tree = shortest_path_tree(res.edge, res.offset);
        return res;
    }

    vector<int> shortest_path_tree(int id, long long offset) const{
        auto [u, v, w] = edges[id];
        vector<long long> from(n, INF);
        from[u] = offset;
        from[v] = 2 * w - offset;

        vector<int> parent_edge = dijkstra(from, 2).second, res;
        for (int x = 0; x < n; x++){
            if (parent_edge[x] != -1) res.push_back(parent_edge[x]);
        }
        /// u and v both stay roots unless one is nearer through the other side, then the center edge joins the halves
        if ((int)res.size() < n - 1) res.push_back(id);

        return res;
    }

    pair<vector<long long>, vector<int>> dijkstra(vector<long long> dist, long long scale) const{
        vector<int> parent_edge(n, -1);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        for (int v = 0; v < n; v++){
            if (dist[v] < INF) heap.push({dist[v], v});
        }

        while (!heap.empty()){
            auto [d, x] = heap.top();
            heap.pop();
            if (d != dist[x]) continue;
            for (auto [y, id] : adj[x]){
                long long nd = d + scale * edges[id].w;
                if (nd < dist[y]) dist[y] = nd, parent_edge[y] = id, heap.push({nd, y});
            }
        }

        return {dist, parent_edge};
    }
};

int main(){
    /***
     *   0 --4-- 1 --6-- 2
     *           |
     *           2
     *           |
     *           3
     * The center sits on edge 1 at distance 1 from vertex 1, so offset 2: radius 5, to vertices 0 and 2
    ***/
    AbsoluteCenter path(4);
    path.add_edge(0, 1, 4), path.add_edge(1, 2, 6), path.add_edge(1, 3, 2);
    AbsoluteCenterResult r = path.solve();
    assert(r.diameter == 10);
    assert(r.edge == 1);
    assert(r.offset == 2);

    /// A triangle of 1, 1, 3: the best spanning tree is the path 0 - 1 - 2, diameter 2, centered on vertex 1
    AbsoluteCenter triangle(3);
    triangle.add_edge(0, 2, 3), triangle.add_edge(0, 1, 1), triangle.add_edge(1, 2, 1);
    AbsoluteCenterResult t = triangle.solve();
    assert(t.diameter == 2);
    sort(t.tree.begin(), t.tree.end());
    assert((t.tree == vector<int>{1, 2}));

    /// A cycle of four unit edges: the center is the middle of any edge, radius 3 / 2 beats 2 at a vertex
    AbsoluteCenter cycle(4);
    for (int i = 0; i < 4; i++) cycle.add_edge(i, (i + 1) % 4, 1);
    AbsoluteCenterResult c = cycle.solve();
    assert(c.diameter == 3);
    assert(c.offset == 1);
    assert(c.tree.size() == 3);

    return 0;
}
