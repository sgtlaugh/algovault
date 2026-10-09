#include "../common.h"

#define main library_main
#include "../../code_library/graphs/complement_graph_bfs.cpp"
#undef main

vector<int> matrix_bfs(const vector<vector<char>>& linked, int s){
    int n = linked.size();
    vector<int> dist(n, -1), frontier = {s};
    dist[s] = 0;
    for (size_t i = 0; i < frontier.size(); i++){
        int u = frontier[i];
        for (int v = 0; v < n; v++){
            if (linked[u][v] && dist[v] == -1) dist[v] = dist[u] + 1, frontier.push_back(v);
        }
    }
    return dist;
}

/// Distances from every source and the component labels against BFS on the explicitly built complement matrix
void check(int n, const vector<pair<int, int>>& edges){
    ComplementGraph g(n);
    vector<vector<char>> linked(n, vector<char>(n, 1));
    for (auto [u, v] : edges) g.add_edge(u, v), linked[u][v] = linked[v][u] = 0;
    for (int v = 0; v < n; v++) linked[v][v] = 0;

    vector<int> comp(n, -1);
    int count = 0;
    for (int s = 0; s < n; s++){
        vector<int> dist = matrix_bfs(linked, s);
        assert(g.bfs(s) == dist);
        if (comp[s] != -1) continue;
        for (int v = 0; v < n; v++){
            if (dist[v] != -1) comp[v] = count;
        }
        count++;
    }
    assert(g.components() == comp);
}

int main(){
    for (int n = 0; n <= 5; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }
        for (int mask = 0; mask < (1 << pairs.size()); mask++){
            vector<pair<int, int>> edges;
            for (size_t i = 0; i < pairs.size(); i++){
                if (mask >> i & 1) edges.push_back(pairs[i]);
            }
            check(n, edges);
        }
    }

    /// Dense originals make sparse complements with long shortest paths, so the edge count spans the whole range
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 60), pairs = n * (n - 1) / 2;
        int extra = stress::rand_int(0, pairs / 4);
        double keep = stress::rand_int(0, 100) / 100.0;
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (stress::rand_int(0, 999) < keep * 1000) edges.push_back({u, v});
            }
        }
        /// Extra random pairs bring in self loops and duplicate edges
        for (int i = 0; i < extra; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            edges.push_back({u, v});
        }
        shuffle(edges.begin(), edges.end(), stress::rng());
        check(n, edges);
    }

    /// Complete graph on 2000 vertices minus a shuffled Hamiltonian path: the complement is that path, about 2e6 original edges
    int n = 2000;
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());
    vector<int> position(n);
    for (int i = 0; i < n; i++) position[order[i]] = i;

    ComplementGraph dense(n);
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (abs(position[u] - position[v]) != 1) dense.add_edge(u, v);
        }
    }
    vector<int> dist = dense.bfs(order[0]);
    for (int v = 0; v < n; v++) assert(dist[v] == position[v]);
    assert((dense.components() == vector<int>(n, 0)));

    /// Sparse original on 3e5 vertices: a star leaves its center isolated in the complement, the rest is one component at distance <= 2
    int big = 300000;
    ComplementGraph star(big);
    for (int v = 1; v < big; v++) star.add_edge(0, v);
    for (int i = 0; i < big; i++) star.add_edge(stress::rand_int(1, big - 1), stress::rand_int(1, big - 1));
    vector<int> comp = star.components();
    assert(comp[0] == 0);
    for (int v = 1; v < big; v++) assert(comp[v] == 1);
    dist = star.bfs(0);
    assert(dist[0] == 0 && count(dist.begin(), dist.end(), -1) == big - 1);
    dist = star.bfs(1);
    assert(dist[0] == -1 && *max_element(dist.begin(), dist.end()) <= 2);
    return 0;
}
