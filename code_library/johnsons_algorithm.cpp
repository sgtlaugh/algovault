/***
 *
 * Johnson's Algorithm
 * All pairs shortest paths on a sparse directed graph with negative edges
 *
 * Complexity: O(n m) for Bellman Ford once, then O(n^2 + n m log n) for n runs of Dijkstra, O(n^2) memory for dist
 *
 * johnson(n, edges, dist): edges[i] = {u, v, w}, 0-based, long long weights (negative allowed)
 *     returns false if the graph has a negative cycle (dist is left empty)
 *     otherwise dist[u][v] is the shortest distance, or JOHNSON_INF if v is unreachable from u
 *     every |path weight| must stay below 1e18: reweighted distances reach up to twice that
 *
 * Bellman Ford from a virtual source gives potentials h, every edge is reweighted to w + h[u] - h[v] >= 0,
 * so plain Dijkstra works from each source and the true distance is recovered as d - h[u] + h[v]
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const long long JOHNSON_INF = numeric_limits<long long>::max() / 4;

bool johnson(int n, const vector<array<long long, 3>>& edges, vector<vector<long long>>& dist){
    dist.clear();
    vector<long long> h(n, 0);  /// the virtual source reaches every vertex with weight 0
    for (int round = 0; round < n; round++){
        bool changed = false;
        for (auto [u, v, w] : edges){
            if (h[u] + w < h[v]) h[v] = h[u] + w, changed = true;
        }
        if (!changed) break;
        if (round == n - 1) return false;  /// shortest paths need at most n - 1 passes, a change in pass n means a negative cycle
    }

    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges) adj[u].push_back({(int)v, w + h[u] - h[v]});

    dist.assign(n, vector<long long>(n, JOHNSON_INF));
    for (int s = 0; s < n; s++){
        vector<long long>& d = dist[s];
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
        d[s] = 0, heap.push({0, s});
        while (!heap.empty()){
            auto [du, u] = heap.top();
            heap.pop();
            if (du != d[u]) continue;
            for (auto [v, w] : adj[u]){
                if (du + w < d[v]) d[v] = du + w, heap.push({d[v], v});
            }
        }
        for (int v = 0; v < n; v++){
            if (d[v] != JOHNSON_INF) d[v] += h[v] - h[s];
        }
    }
    return true;
}

int main(){
    vector<vector<long long>> dist;
    assert(johnson(3, {{0, 1, 2}, {1, 2, -15}, {0, 2, -10}}, dist));
    assert(dist[0][1] == 2 && dist[0][2] == -13 && dist[1][2] == -15 && dist[0][0] == 0);
    assert(dist[2][0] == JOHNSON_INF && dist[1][0] == JOHNSON_INF);

    assert(!johnson(3, {{0, 1, 1}, {1, 2, -3}, {2, 0, 1}}, dist) && dist.empty());
    assert(!johnson(1, {{0, 0, -1}}, dist));
    assert(johnson(1, {{0, 0, 5}}, dist) && dist[0][0] == 0);

    assert(johnson(4, {{0, 1, 5}, {0, 1, 3}, {1, 2, -2}, {2, 3, 4}, {3, 1, -2}}, dist));
    assert(dist[0][3] == 5 && dist[3][2] == -4 && dist[2][1] == 2);
    return 0;
}
