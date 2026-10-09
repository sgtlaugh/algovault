#include "../common.h"

#define main library_main
#include "../../code_library/graphs/range_edge_dijkstra.cpp"
#undef main

const long long INF = RangeEdgeGraph::INF;

/// O(n^2) Dijkstra on an adjacency matrix holding every range edge expanded into single edges
vector<long long> expanded_dijkstra(const vector<vector<long long>>& w, int s){
    int n = w.size();
    vector<long long> dist(n, LLONG_MAX);
    vector<bool> done(n, false);
    dist[s] = 0;

    for (int round = 0; round < n; round++){
        int u = -1;
        for (int v = 0; v < n; v++){
            if (!done[v] && dist[v] != LLONG_MAX && (u == -1 || dist[v] < dist[u])) u = v;
        }
        if (u == -1 || dist[u] >= INF) break;  /// stops before dist[u] + w can overflow
        done[u] = true;
        for (int v = 0; v < n; v++){
            if (w[u][v] != LLONG_MAX) dist[v] = min(dist[v], dist[u] + w[u][v]);
        }
    }

    for (auto& d : dist) d = min(d, INF);
    return dist;
}

long long random_weight(int mode){
    if (mode == 0) return stress::rand_int(0, 5);
    if (mode == 1) return stress::rand_int(0, 1000000000000LL);
    return stress::rand_int(INF / 3, INF);
}

void add_random_edges(RangeEdgeGraph& g, vector<vector<long long>>& w, int count, int mode){
    int n = w.size();
    for (int i = 0; i < count; i++){
        int type = stress::rand_int(0, 2), u = stress::rand_int(0, n - 1);
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        long long c = random_weight(mode);

        if (type == 0){
            g.add_edge(u, l, c);
            w[u][l] = min(w[u][l], c);
        }
        else if (type == 1){
            g.add_edge_to_range(u, l, r, c);
            for (int v = l; v <= r; v++) w[u][v] = min(w[u][v], c);
        }
        else{
            g.add_edge_from_range(l, r, u, c);
            for (int v = l; v <= r; v++) w[v][u] = min(w[v][u], c);
        }
    }
}

/// Every i reaches the next k vertices through one range edge, so vertex v is ceil(v / k) hops away
void check_long_chains(int n, int k){
    RangeEdgeGraph to_range(n), from_range(n);
    for (int i = 0; i < n; i++){
        if (i + 1 < n) to_range.add_edge_to_range(i, i + 1, min(n - 1, i + k), 1);
        if (i > 0) from_range.add_edge_from_range(max(0, i - k), i - 1, i, 1);
    }

    vector<long long> a = to_range.dijkstra(0), b = from_range.dijkstra(0);
    for (int v = 0; v < n; v++) assert(a[v] == (v + k - 1) / k && b[v] == (v + k - 1) / k);
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 600 ? it % 12 + 1 : stress::rand_int(1, 70);
        int mode = it % 10 == 9 ? 2 : it % 3 == 2;
        RangeEdgeGraph g(n);
        vector<vector<long long>> w(n, vector<long long>(n, LLONG_MAX));

        for (int batch = 0; batch < 3; batch++){
            add_random_edges(g, w, stress::rand_int(0, 2 * n), mode);
            int s = stress::rand_int(0, n - 1);
            assert(g.dijkstra(s) == expanded_dijkstra(w, s));
        }
    }

    check_long_chains(100000, 1);
    check_long_chains(100000, 777);
    check_long_chains(99999, 99999);

    /// Distances reaching INF exactly are unreachable, one below INF is kept
    RangeEdgeGraph edge_of_infinity(4);
    edge_of_infinity.add_edge_to_range(0, 1, 3, INF / 2);
    edge_of_infinity.add_edge_from_range(1, 1, 2, INF / 2);
    edge_of_infinity.add_edge_from_range(1, 1, 3, INF / 2 - 1);
    edge_of_infinity.add_edge(3, 0, INF / 2);
    edge_of_infinity.add_edge(2, 0, INF);
    assert((edge_of_infinity.dijkstra(0) == vector<long long>{0, INF / 2, INF / 2, INF / 2}));
    assert((edge_of_infinity.dijkstra(1) == vector<long long>{INF - 1, 0, INF / 2, INF / 2 - 1}));
    assert((edge_of_infinity.dijkstra(2) == vector<long long>{INF, INF, 0, INF}));
    assert((edge_of_infinity.dijkstra(3) == vector<long long>{INF / 2, INF, INF, 0}));
    return 0;
}
