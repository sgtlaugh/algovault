/***
 *
 * Dijkstra's Algorithm and 0-1 BFS
 * Shortest paths from one or more sources with non-negative edge weights, plus shortest path reconstruction
 *
 * Complexity:
 *   - Dijkstra: O((n + m) log m) per run with a binary heap
 *   - DenseDijkstra: O(n^2) per run, O(n^2) memory for the weight matrix (32 MB at n = 2000)
 *   - ZeroOneBFS: O(n + m) per run with a deque
 *
 * Dijkstra, DenseDijkstra and ZeroOneBFS share one API on 0-based vertices, shown for Dijkstra:
 * Dijkstra g(n): an empty graph on n vertices
 * g.add_edge(u, v, w): adds the directed edge u -> v of weight w, pass directed = false for an undirected one
 * g.run(src): fills g.dist and g.parent with the shortest paths from src
 * g.run(src, target): stops as soon as target is settled, so only dist[target] and path(target) are final
 * g.run(sources): the same from the nearest of several sources, all at distance 0, and also takes a target
 * g.path(v): the vertices of one shortest path from a source to v, empty if v is unreachable
 * g.dist[v]: the shortest distance from a source to v, or DIJKSTRA_INF if v is unreachable
 *
 * Notes:
 *   - Weights are long long in [0, DIJKSTRA_INF), 0 or 1 for ZeroOneBFS
 *   - Every shortest distance must stay below DIJKSTRA_INF (~2.3e18)
 *   - Parallel edges and self loops are allowed, and add_edge and run can be called in any order
 *   - DenseDijkstra beats Dijkstra once m passes about n^2 / 2, or when most edges relax, like on a complete DAG
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long DIJKSTRA_INF = numeric_limits<long long>::max() / 4;

vector<int> trace_path(const vector<long long>& dist, const vector<int>& parent, int v){
    if (dist[v] == DIJKSTRA_INF) return {};

    vector<int> path;
    for (; v != -1; v = parent[v]) path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

struct Dijkstra{
    int n;
    vector<tuple<int, int, long long>> arcs;
    vector<int> start;
    vector<pair<int, long long>> adj;  /// u's arcs are adj[start[u] .. start[u + 1]), rebuilt by run after add_edge
    vector<long long> dist;
    vector<int> parent;

    Dijkstra(int n) : n(n), start(n + 1), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w >= 0);
        arcs.push_back({u, v, w});
        if (!directed) arcs.push_back({v, u, w});
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        if (adj.size() != arcs.size()) build_adjacency();
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;

        for (int src : sources) dist[src] = 0, heap.push({0, src});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != dist[u]) continue;
            if (u == target) break;
            for (auto it = adj.begin() + start[u], end = adj.begin() + start[u + 1]; it != end; it++){
                auto [v, w] = *it;
                if (du + w < dist[v]) dist[v] = du + w, parent[v] = u, heap.push({dist[v], v});
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }

    void build_adjacency(){
        fill(start.begin(), start.end(), 0);
        for (auto [u, v, w] : arcs) start[u + 1]++;
        for (int u = 0; u < n; u++) start[u + 1] += start[u];

        adj.resize(arcs.size());
        vector<int> next_slot(start.begin(), start.end() - 1);
        for (auto [u, v, w] : arcs) adj[next_slot[u]++] = {v, w};
    }
};

struct DenseDijkstra{
    int n;
    vector<long long> weight;  /// weight[u * n + v], the lightest u -> v edge or DIJKSTRA_INF
    vector<long long> dist;
    vector<int> parent;

    DenseDijkstra(int n) : n(n), weight((size_t)n * n, DIJKSTRA_INF), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w >= 0);
        long long& cell = weight[(size_t)u * n + v];
        cell = min(cell, w);
        if (!directed) add_edge(v, u, w, true);
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        vector<bool> done(n, false);

        for (int src : sources) dist[src] = 0;
        for (int it = 0; it < n; it++){
            int u = -1;
            for (int v = 0; v < n; v++){
                if (!done[v] && (u == -1 || dist[v] < dist[u])) u = v;
            }
            if (dist[u] == DIJKSTRA_INF || u == target) break;

            done[u] = true;
            /// A missing edge holds DIJKSTRA_INF, so dist[u] + INF >= INF >= dist[v] never relaxes and needs no branch
            const long long* row = &weight[(size_t)u * n];
            for (int v = 0; v < n; v++){
                if (dist[u] + row[v] < dist[v]) dist[v] = dist[u] + row[v], parent[v] = u;
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }
};

struct ZeroOneBFS{
    int n;
    vector<vector<pair<int, int>>> adj;
    vector<long long> dist;
    vector<int> parent;

    ZeroOneBFS(int n) : n(n), adj(n), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w == 0 || w == 1);
        adj[u].push_back({v, (int)w});
        if (!directed) adj[v].push_back({u, (int)w});
    }

    void run(int src, int target = -1){
        run(vector<int>{src}, target);
    }

    void run(const vector<int>& sources, int target = -1){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        deque<int> dq;

        /// The deque only ever holds distances d and d + 1, so a vertex improves at most once after it is first reached
        for (int src : sources) dist[src] = 0, dq.push_back(src);
        while (!dq.empty()){
            int u = dq.front();
            dq.pop_front();
            if (u == target) break;
            for (auto [v, w] : adj[u]){
                if (dist[u] + w >= dist[v]) continue;
                dist[v] = dist[u] + w, parent[v] = u;
                if (w == 0) dq.push_front(v);
                else dq.push_back(v);
            }
        }
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }
};

int main(){
    /// 0 -> 2 -> 1 costs 3 and beats the direct edge of 4, and nothing reaches 4
    Dijkstra g(5);
    g.add_edge(0, 1, 4);
    g.add_edge(0, 2, 1);
    g.add_edge(2, 1, 2);
    g.add_edge(1, 3, 5, false);
    g.run(0);
    assert((g.dist == vector<long long>{0, 3, 1, 8, DIJKSTRA_INF}));
    assert((g.path(3) == vector<int>{0, 2, 1, 3}));
    assert(g.path(4).empty());

    g.run(3, 1);
    assert(g.dist[1] == 5);
    assert((g.path(1) == vector<int>{3, 1}));

    g.run({2, 3});
    assert((g.dist == vector<long long>{DIJKSTRA_INF, 2, 0, 0, DIJKSTRA_INF}));
    assert((g.path(1) == vector<int>{2, 1}));

    /// The lighter of two parallel edges wins
    DenseDijkstra dense(3);
    dense.add_edge(0, 1, 7);
    dense.add_edge(0, 1, 2);
    dense.add_edge(1, 2, 3, false);
    dense.run(0);
    assert((dense.dist == vector<long long>{0, 2, 5}));

    /// 0 -> 2 -> 3 costs 1 through the 0 edge, against 2 through vertex 1
    ZeroOneBFS bfs(4);
    bfs.add_edge(0, 1, 1);
    bfs.add_edge(1, 3, 1);
    bfs.add_edge(0, 2, 1);
    bfs.add_edge(2, 3, 0);
    bfs.run(0);
    assert((bfs.dist == vector<long long>{0, 1, 1, 1}));
    assert((bfs.path(3) == vector<int>{0, 2, 3}));

    return 0;
}
