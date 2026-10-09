/***
 *
 * Dijkstra's Algorithm and 0-1 BFS
 * Single source shortest paths with non-negative edge weights, plus shortest path reconstruction
 *
 * Complexity:
 *     Dijkstra       O((n + m) log m) with a binary heap, O(n + m) memory
 *     DenseDijkstra  O(n^2) per run, O(n^2) memory for the weight matrix (32 MB at n = 2000)
 *     ZeroOneBFS     O(n + m) with a deque, weights must be 0 or 1
 *
 * All three structs share the same API, 0-based vertices:
 *     Graph g(n);
 *     g.add_edge(u, v, w);           directed edge u -> v, pass directed = false for an undirected one
 *     g.run(src);                    fills g.dist and g.parent
 *     g.dist[v]                      shortest distance from src, or DIJKSTRA_INF if v is unreachable
 *     g.path(v)                      vertices of one shortest path src, ..., v, empty if v is unreachable
 *
 * Weights are long long in [0, DIJKSTRA_INF) and every shortest distance must stay below DIJKSTRA_INF (~2.3e18)
 * Parallel edges and self loops are allowed, run can be called again with another source on the same graph
 *
 * DenseDijkstra is for dense graphs where the heap version hits its worst case, measured at -O2 with n = 2000:
 *     complete DAG weighted so every edge relaxes (2e6 heap pushes): DenseDijkstra 17 ms, Dijkstra 685 ms
 *     random weights: Dijkstra 12x faster with 2% of all edges present, down to 1.3x with 100%, since few edges relax
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
    vector<vector<pair<int, long long>>> adj;
    vector<long long> dist;
    vector<int> parent;

    Dijkstra(int n) : n(n), adj(n), dist(n), parent(n) {}

    void add_edge(int u, int v, long long w, bool directed = true){
        assert(w >= 0);
        adj[u].push_back({v, w});
        if (!directed) adj[v].push_back({u, w});
    }

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }

    void run(int src){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;

        dist[src] = 0, heap.push({0, src});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != dist[u]) continue;  /// stale entry, u was already settled with a smaller distance
            for (auto [v, w] : adj[u]){
                if (du + w < dist[v]) dist[v] = du + w, parent[v] = u, heap.push({dist[v], v});
            }
        }
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

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }

    void run(int src){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        vector<bool> done(n, false);

        dist[src] = 0;
        for (int it = 0; it < n; it++){
            int u = -1;
            for (int v = 0; v < n; v++){
                if (!done[v] && (u == -1 || dist[v] < dist[u])) u = v;
            }
            if (dist[u] == DIJKSTRA_INF) break;

            done[u] = true;
            /// A missing edge holds DIJKSTRA_INF, so dist[u] + INF >= INF >= dist[v] never relaxes and needs no branch
            const long long* row = &weight[(size_t)u * n];
            for (int v = 0; v < n; v++){
                if (dist[u] + row[v] < dist[v]) dist[v] = dist[u] + row[v], parent[v] = u;
            }
        }
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

    vector<int> path(int v) const{
        return trace_path(dist, parent, v);
    }

    void run(int src){
        fill(dist.begin(), dist.end(), DIJKSTRA_INF);
        fill(parent.begin(), parent.end(), -1);
        deque<int> dq;

        /// The deque only ever holds distances d and d + 1, so a vertex improves at most once after it is first reached
        dist[src] = 0, dq.push_back(src);
        while (!dq.empty()){
            int u = dq.front();
            dq.pop_front();
            for (auto [v, w] : adj[u]){
                if (dist[u] + w >= dist[v]) continue;
                dist[v] = dist[u] + w, parent[v] = u;
                if (w == 0) dq.push_front(v);
                else dq.push_back(v);
            }
        }
    }
};

template<typename Graph>
void self_test(){
    Graph g(6);
    g.add_edge(0, 1, 1);
    g.add_edge(0, 2, 1);
    g.add_edge(1, 3, 1, false);
    g.add_edge(2, 3, 0);
    g.add_edge(3, 4, 1);
    g.add_edge(4, 4, 0);
    g.run(0);
    assert((g.dist == vector<long long>{0, 1, 1, 1, 2, DIJKSTRA_INF}));
    assert((g.path(4) == vector<int>{0, 2, 3, 4}));
    assert((g.path(0) == vector<int>{0}));
    assert(g.path(5).empty());

    g.run(3);
    assert((g.dist == vector<long long>{DIJKSTRA_INF, 1, DIJKSTRA_INF, 0, 1, DIJKSTRA_INF}));
    assert((g.path(1) == vector<int>{3, 1}));

    Graph single(1);
    single.run(0);
    assert((single.dist == vector<long long>{0}) && (single.path(0) == vector<int>{0}));
}

template<typename Graph>
void weighted_test(){
    Graph g(5);
    g.add_edge(0, 1, 10);
    g.add_edge(0, 1, 4);
    g.add_edge(0, 2, 1);
    g.add_edge(2, 1, 2);
    g.add_edge(1, 3, 5, false);
    g.add_edge(2, 3, 8);
    g.add_edge(3, 4, 3);
    g.run(0);
    assert((g.dist == vector<long long>{0, 3, 1, 8, 11}));
    assert((g.path(4) == vector<int>{0, 2, 1, 3, 4}));

    g.run(3);
    assert((g.dist == vector<long long>{DIJKSTRA_INF, 5, DIJKSTRA_INF, 0, 3}));

    const long long big = 1000000000000000000LL;
    Graph far(3);
    far.add_edge(0, 1, big);
    far.add_edge(1, 2, big);
    far.add_edge(0, 2, DIJKSTRA_INF - 1);
    far.run(0);
    assert((far.dist == vector<long long>{0, big, 2 * big}));
    assert((far.path(2) == vector<int>{0, 1, 2}));
}

int main(){
    self_test<Dijkstra>();
    self_test<DenseDijkstra>();
    self_test<ZeroOneBFS>();

    weighted_test<Dijkstra>();
    weighted_test<DenseDijkstra>();

    return 0;
}
