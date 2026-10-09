#include "../common.h"

#define main library_main
#include "../../code_library/graphs/small_cycle_counting.cpp"
#undef main

/// O(n^3) triangles over vertex triples, O(n^4) 4-cycles over ordered 4-tuples where each cycle appears 8 times
pair<long long, long long> brute(int n, const vector<vector<char>>& has){
    long long tri = 0, quad = 0;

    for (int a = 0; a < n; a++){
        for (int b = a + 1; b < n; b++){
            for (int c = b + 1; c < n; c++) tri += has[a][b] && has[b][c] && has[a][c];
        }
    }

    for (int a = 0; a < n; a++){
        for (int b = 0; b < n; b++){
            if (!has[a][b]) continue;
            for (int c = 0; c < n; c++){
                if (c == a || !has[b][c]) continue;
                for (int d = 0; d < n; d++) quad += d != b && has[c][d] && has[d][a];
            }
        }
    }

    return {tri, quad / 8};
}

void check(int n, const vector<pair<int, int>>& edges){
    vector<vector<char>> has(n, vector<char>(n, 0));
    SmallCycles g(n);
    for (auto [u, v] : edges){
        has[u][v] = has[v][u] = 1;
        g.add_edge(u, v);
    }

    auto [tri, quad] = brute(n, has);
    assert(g.triangles() == tri && g.four_cycles() == quad);
    assert(g.triangles() == tri && g.four_cycles() == quad);
}

/// Random density plus optional hubs joined to most vertices, so high degree vertices get every index
vector<pair<int, int>> random_graph(int n){
    vector<vector<char>> has(n, vector<char>(n, 0));
    int density = stress::rand_int(0, 100), hubs = stress::rand_int(0, 3), hub_density = stress::rand_int(50, 100);

    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            bool hub = u < hubs || v < hubs;
            has[u][v] = stress::rand_int(0, 99) < (hub ? hub_density : density);
        }
    }

    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (!has[u][v]) continue;
            edges.push_back({label[u], label[v]});
            if (stress::rand_int(0, 1)) swap(edges.back().first, edges.back().second);
        }
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

int main(){
    /// Every labeled simple graph on up to 6 vertices
    for (int n = 0; n <= 6; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }
        for (int mask = 0; mask < (1 << pairs.size()); mask++){
            vector<pair<int, int>> edges;
            for (int i = 0; i < (int)pairs.size(); i++){
                if (mask >> i & 1) edges.push_back(pairs[i]);
            }
            check(n, edges);
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 14 : 40);
        check(n, random_graph(n));
    }

    /// Two hubs of degree 1e5 on the lowest indices: ordering by index alone scans 1e10 paths
    int leaves = 100000;
    SmallCycles book(leaves + 2);
    for (int i = 2; i < leaves + 2; i++) book.add_edge(0, i), book.add_edge(1, i);

    auto start = chrono::steady_clock::now();
    assert(book.triangles() == 0 && book.four_cycles() == 1LL * leaves * (leaves - 1) / 2);
    assert(chrono::steady_clock::now() - start < chrono::seconds(5));

    /// K_400, m = 79800: C(n, 3) triangles and 3 C(n, 4) 4-cycles
    long long n = 400;
    SmallCycles dense(n);
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++) dense.add_edge(u, v);
    }

    start = chrono::steady_clock::now();
    assert(dense.triangles() == n * (n - 1) * (n - 2) / 6);
    assert(dense.four_cycles() == n * (n - 1) * (n - 2) * (n - 3) / 8);
    assert(chrono::steady_clock::now() - start < chrono::seconds(5));
    return 0;
}
