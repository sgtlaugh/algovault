/***
 *
 * Minimum Weight Cycle
 * The lightest simple cycle of a directed or undirected graph with non-negative weights, unweighted girth included
 *
 * Complexity:
 *   - O(n + n m log m) with a Dijkstra from every vertex
 *   - O(n + n m) with a BFS from every vertex when all weights are equal, the girth case
 *
 * MinimumWeightCycle g(n, MinimumWeightCycle::DIRECTED): empty graph on vertices 0 to n - 1, or UNDIRECTED
 * g.add_edge(u, v, w): edge u -> v, or u - v, of weight w >= 0, w = 1 by default
 * g.solve(): {exists, weight, vertices, edges}, a lightest simple cycle
 *
 * Notes:
 *   - vertices and edges list the cycle in traversal order, edges as ids numbered from 0 in add_edge order
 *   - vertices[i] -> vertices[i + 1] uses edges[i], the last edge returns to vertices[0]
 *   - exists is false when the graph has no cycle, then weight = 0 and both lists are empty
 *   - A self loop is a cycle of one; undirected cycles never reuse an edge, so two parallel edges are a cycle of two
 *   - Cycles of 3 or more vertices only: keep the lightest of each set of parallel undirected edges, drop self loops
 *   - The weight of every path must stay below 1e18, so cycles weigh below 2e18
 *   - Each search stops at the lightest cycle found so far, half of it when undirected: fast when short cycles exist
 *   - Dense directed graph: floyd_warshall.cpp, then the min over edges v -> u of dist[u][v] + w, is O(n^3)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct WeightCycle{
    bool exists;
    long long weight;
    vector<int> vertices, edges;
};

struct MinimumWeightCycle{
    enum Kind{DIRECTED, UNDIRECTED};
    static constexpr long long INF = 4e18;  /// above two paths of 1e18 plus an edge, so every cycle is lighter

    struct Edge{
        int u, v;
        long long w;
    };

    int n;
    Kind kind;
    vector<Edge> edges;
    vector<vector<pair<int, int>>> adj;
    vector<long long> dist;
    vector<int> parent_edge, top, touched;
    vector<char> done;

    MinimumWeightCycle(int n, Kind kind) : n(n), kind(kind), adj(n) {}

    void add_edge(int u, int v, long long w = 1){
        assert(w >= 0);
        adj[u].push_back({v, (int)edges.size()});
        if (kind == UNDIRECTED) adj[v].push_back({u, (int)edges.size()});
        edges.push_back({u, v, w});
    }

    WeightCycle solve(){
        bool equal_weights = true;
        for (Edge& e : edges) equal_weights &= e.w == edges[0].w;

        long long best = INF;
        int best_source = -1;
        for (int s = 0; s < n; s++){
            long long weight = search(s, equal_weights, kind == DIRECTED ? best : (best + 1) / 2).first;
            if (weight < best) best = weight, best_source = s;
        }
        if (best_source == -1) return {false, 0, {}, {}};

        int best_edge = search(best_source, equal_weights, INF).second;
        WeightCycle res{true, best, {}, {}};
        auto [u, v, w] = edges[best_edge];
        for (int x = u; x != best_source; x = other_end(parent_edge[x], x)){
            res.vertices.push_back(x);
            res.edges.push_back(parent_edge[x]);
        }
        res.vertices.push_back(best_source);
        reverse(res.vertices.begin(), res.vertices.end());
        reverse(res.edges.begin(), res.edges.end());

        res.edges.push_back(best_edge);
        for (int x = v; x != best_source; x = other_end(parent_edge[x], x)){
            res.vertices.push_back(x);
            res.edges.push_back(parent_edge[x]);
        }

        return res;
    }

    /// The lightest cycle {weight, closing edge} through s found settling only dist < limit, INF if none
    pair<long long, int> search(int s, bool equal_weights, long long limit){
        /// Resetting only what the last search reached makes a search cost what it visits, not O(n) on sparse graphs
        dist.resize(n, INF), parent_edge.resize(n, -1), top.resize(n, -1), done.resize(n, 0);
        for (int v : touched) dist[v] = INF, parent_edge[v] = -1, top[v] = -1, done[v] = 0;
        touched.assign(1, s);
        dist[s] = 0, top[s] = s;
        pair<long long, int> res{INF, -1};

        /// With equal weights, first-in first-out order is distance order, so a queue replaces the heap
        vector<int> queue;
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        if (equal_weights) queue.push_back(s);
        else heap.push({0, s});
        for (size_t head = 0; head < queue.size() || !heap.empty();){
            int u;
            if (equal_weights) u = queue[head++];
            else{
                auto [d, x] = heap.top();
                heap.pop();
                if (d != dist[x]) continue;
                u = x;
            }
            if (dist[u] >= limit) break;  /// every vertex of a cycle through s lighter than best is closer than limit
            done[u] = 1;

            for (auto [v, id] : adj[u]){
                long long w = edges[id].w;
                if (done[v]){
                    if (kind == DIRECTED && v != s) continue;

                    /// top is the child of s above each vertex: paths into different subtrees meet only at s
                    if (kind == UNDIRECTED && u != v && (parent_edge[u] == id || top[u] == top[v])) continue;
                    if (kind == UNDIRECTED && u == v && u != s) continue;

                    res = min(res, {dist[u] + dist[v] + w, id});
                    continue;
                }

                if (dist[u] + w >= dist[v]) continue;
                if (dist[v] == INF) touched.push_back(v);
                dist[v] = dist[u] + w;
                parent_edge[v] = id;
                top[v] = u == s ? v : top[u];
                if (equal_weights) queue.push_back(v);
                else heap.push({dist[v], v});
            }
        }

        return res;
    }

    int other_end(int id, int x) const{
        return edges[id].u == x ? edges[id].v : edges[id].u;
    }
};

int main(){
    /// Undirected: triangle 0-1-2 of weight 3 + 4 + 5 = 12 and a square 2-3-4-5 of weight 1 + 1 + 1 + 1 = 4
    MinimumWeightCycle g(6, MinimumWeightCycle::UNDIRECTED);
    g.add_edge(0, 1, 3);
    g.add_edge(1, 2, 4);
    g.add_edge(2, 0, 5);
    g.add_edge(2, 3, 1);
    g.add_edge(3, 4, 1);
    g.add_edge(4, 5, 1);
    g.add_edge(5, 2, 1);
    WeightCycle c = g.solve();
    assert(c.exists);
    assert(c.weight == 4);
    assert(c.vertices.size() == 4);

    /// Two parallel edges 0 - 1 of weight 2 and 3 are a cycle of weight 5, one edge alone is not a cycle
    MinimumWeightCycle pair_of(2, MinimumWeightCycle::UNDIRECTED);
    pair_of.add_edge(0, 1, 2);
    assert(!pair_of.solve().exists);
    pair_of.add_edge(1, 0, 3);
    WeightCycle two = pair_of.solve();
    assert(two.weight == 5);
    assert((two.vertices == vector<int>{0, 1}));
    assert((two.edges == vector<int>{0, 1}));

    /// Directed: 0 -> 1 -> 2 -> 0 weighs 6, beating the 2-cycle 1 -> 2 -> 1 of weight 2 + 5
    MinimumWeightCycle d(3, MinimumWeightCycle::DIRECTED);
    d.add_edge(0, 1, 1);
    d.add_edge(1, 2, 2);
    d.add_edge(2, 0, 3);
    d.add_edge(2, 1, 5);
    WeightCycle dc = d.solve();
    assert(dc.weight == 6);
    assert((dc.vertices == vector<int>{0, 1, 2}));
    assert((dc.edges == vector<int>{0, 1, 2}));

    /// Unweighted girth: the 5-cycle 0-1-2-3-4 with a chord 0-2 has girth 3
    MinimumWeightCycle pentagon(5, MinimumWeightCycle::UNDIRECTED);
    for (int i = 0; i < 5; i++) pentagon.add_edge(i, (i + 1) % 5);
    assert(pentagon.solve().weight == 5);
    pentagon.add_edge(0, 2);
    assert(pentagon.solve().weight == 3);

    return 0;
}
