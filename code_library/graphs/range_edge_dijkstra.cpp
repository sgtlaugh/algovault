/***
 *
 * Dijkstra with Range Edges
 * Single source shortest paths where an edge can go from a vertex to a whole range, or from a whole range to a vertex
 *
 * Complexity: O(n + m log n) graph nodes and edges, O((n + m log n) log n) per dijkstra call, m = number of add calls
 *
 * RangeEdgeGraph g(n): n >= 1, vertices 0 to n - 1, ranges [l, r] are inclusive and require 0 <= l <= r < n
 *     g.add_edge(u, v, w): u -> v
 *     g.add_edge_to_range(u, l, r, w): u -> every v in [l, r]
 *     g.add_edge_from_range(l, r, v, w): every u in [l, r] -> v
 *     g.dijkstra(s) returns the distances of the n vertices from s, RangeEdgeGraph::INF if unreachable
 * Weights must be in [0, INF], INF = 4e18, a vertex whose shortest distance is >= INF is reported as INF
 * Several dijkstra calls can reuse one graph, edges can be added between calls
 *
 * Two segment trees over the vertices share the vertices as leaves: in the out-tree every node has
 * 0 weight edges down to its children, in the in-tree up to its parent. A range edge splits [l, r] into
 * O(log n) tree nodes and links u to out-tree nodes (or in-tree nodes to v). The trees are separate so
 * a range edge never lets one vertex of the range reach another for free (CF 786B Legacy)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct RangeEdgeGraph{
    static constexpr long long INF = 4000000000000000000LL;

    int n;
    vector<vector<pair<int, long long>>> adj;

    /// Tree node i in [1, n) is graph node n + i in the out-tree and 2n + i in the in-tree, leaf n + v is vertex v
    RangeEdgeGraph(int n): n(n), adj(3 * n){
        for (int i = 1; i < n; i++){
            for (int c : {2 * i, 2 * i + 1}){
                adj[out_node(i)].push_back({out_node(c), 0});
                adj[in_node(c)].push_back({in_node(i), 0});
            }
        }
    }

    void add_edge(int u, int v, long long w){
        adj[u].push_back({v, w});
    }

    void add_edge_from_range(int l, int r, int v, long long w){
        for (int i : cover(l, r)) adj[in_node(i)].push_back({v, w});
    }

    void add_edge_to_range(int u, int l, int r, long long w){
        for (int i : cover(l, r)) adj[u].push_back({out_node(i), w});
    }

    vector<long long> dijkstra(int s){
        vector<long long> dist(3 * n, INF);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        dist[s] = 0, heap.push({0, s});

        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != dist[u]) continue;
            for (auto [v, w] : adj[u]){
                if (du + w < dist[v]) dist[v] = du + w, heap.push({dist[v], v});
            }
        }

        dist.resize(n);
        return dist;
    }

private:
    /// Bottom-up segment tree nodes whose leaves are exactly [l, r], valid for any n, not only powers of two
    vector<int> cover(int l, int r) const{
        vector<int> nodes;
        for (l += n, r += n + 1; l < r; l >>= 1, r >>= 1){
            if (l & 1) nodes.push_back(l++);
            if (r & 1) nodes.push_back(--r);
        }
        return nodes;
    }

    int in_node(int i) const{
        return i >= n ? i - n : 2 * n + i;
    }

    int out_node(int i) const{
        return i >= n ? i - n : n + i;
    }
};

int main(){
    const long long INF = RangeEdgeGraph::INF;

    RangeEdgeGraph sample1(3);
    sample1.add_edge_to_range(2, 1, 2, 17);
    sample1.add_edge_to_range(2, 1, 1, 16);
    sample1.add_edge_to_range(1, 1, 2, 3);
    sample1.add_edge_from_range(0, 0, 2, 12);
    sample1.add_edge(2, 2, 17);
    assert((sample1.dijkstra(0) == vector<long long>{0, 28, 12}));

    RangeEdgeGraph sample2(4);
    sample2.add_edge_from_range(0, 2, 3, 12);
    sample2.add_edge_to_range(1, 2, 3, 10);
    sample2.add_edge(1, 3, 16);
    assert((sample2.dijkstra(0) == vector<long long>{0, INF, INF, 12}));
    assert((sample2.dijkstra(1) == vector<long long>{INF, 0, 10, 10}));

    RangeEdgeGraph single(1);
    single.add_edge_to_range(0, 0, 0, 5);
    single.add_edge_from_range(0, 0, 0, 7);
    assert((single.dijkstra(0) == vector<long long>{0}));

    RangeEdgeGraph leak(5);
    leak.add_edge_to_range(0, 1, 3, 5);
    leak.add_edge_from_range(1, 3, 4, 2);
    assert((leak.dijkstra(0) == vector<long long>{0, 5, 5, 5, 7}));
    assert((leak.dijkstra(2) == vector<long long>{INF, INF, 0, INF, 2}));
    assert((leak.dijkstra(4) == vector<long long>{INF, INF, INF, INF, 0}));

    RangeEdgeGraph mixed(6);
    mixed.add_edge_to_range(0, 2, 5, 10);
    mixed.add_edge(0, 1, 1);
    mixed.add_edge_from_range(0, 1, 4, 2);
    mixed.add_edge_to_range(4, 5, 5, 0);
    mixed.add_edge(5, 3, 1);
    assert((mixed.dijkstra(0) == vector<long long>{0, 1, 10, 3, 2, 2}));
    mixed.add_edge(1, 2, 3);
    assert((mixed.dijkstra(0) == vector<long long>{0, 1, 4, 3, 2, 2}));
    return 0;
}
