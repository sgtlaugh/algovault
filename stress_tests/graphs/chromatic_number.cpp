#include "../common.h"

#define main library_main
#include "../../code_library/graphs/chromatic_number.cpp"
#undef main

/// Vertex v may open at most one new color, so color permutations are not searched again
bool can_color(int v, int used, int k, const vector<vector<bool>>& adj, vector<int>& color){
    int n = adj.size();
    if (v == n) return true;

    for (int c = 0; c < min(k, used + 1); c++){
        bool clash = false;
        for (int u = 0; u < v; u++) clash |= adj[v][u] && color[u] == c;
        if (clash) continue;
        color[v] = c;
        if (can_color(v + 1, max(used, c + 1), k, adj, color)) return true;
    }
    return false;
}

/// Smallest k for which backtracking finds a proper coloring
int brute(const vector<vector<bool>>& adj){
    int n = adj.size();
    vector<int> color(n);
    for (int k = 0;; k++){
        if (can_color(0, 0, k, adj, color)) return k;
    }
}

int solve(int n, const vector<pair<int, int>>& edges){
    ChromaticNumber g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);
    return g.solve();
}

/// Random proper c-coloring with every color used, a rainbow c-clique and up to 3n random edges between colors: chromatic number exactly c
int planted(int n, int c){
    vector<int> color(n);
    for (int v = 0; v < n; v++) color[v] = v < c ? v : stress::rand_int(0, c - 1);
    shuffle(color.begin(), color.end(), stress::rng());

    vector<int> first(c, -1);
    for (int v = 0; v < n; v++){
        if (first[color[v]] == -1) first[color[v]] = v;
    }

    vector<pair<int, int>> edges;
    for (int a = 0; a < c; a++){
        for (int b = a + 1; b < c; b++) edges.push_back({first[a], first[b]});
    }
    int attempts = stress::rand_int(0, 3 * n);
    for (int i = 0; i < attempts; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (color[u] != color[v]) edges.push_back({u, v});
    }

    shuffle(edges.begin(), edges.end(), stress::rng());
    return solve(n, edges);
}

int main(){
    /// Every labeled graph on up to 5 vertices
    for (int n = 0; n <= 5; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }

        for (int bits = 0; bits < (1 << pairs.size()); bits++){
            vector<vector<bool>> adj(n, vector<bool>(n));
            vector<pair<int, int>> edges;
            for (int i = 0; i < (int)pairs.size(); i++){
                if (!(bits >> i & 1)) continue;
                auto [u, v] = pairs[i];
                adj[u][v] = adj[v][u] = true;
                edges.push_back(stress::rand_int(0, 1) ? pairs[i] : make_pair(v, u));
            }
            assert(solve(n, edges) == brute(adj));
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(6, 9);
        int density = stress::rand_int(1, 99);
        vector<vector<bool>> adj(n, vector<bool>(n));
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (stress::rand_int(1, 100) > density) continue;
                adj[u][v] = adj[v][u] = true;
                edges.push_back({u, v});
                if (stress::rand_int(0, 9) == 0) edges.push_back({v, u});
            }
        }
        assert(solve(n, edges) == brute(adj));
    }

    for (long long it = 0; it < stress::scaled(40); it++){
        int n = stress::rand_int(10, 16), c = stress::rand_int(1, n);
        assert(planted(n, c) == c);
    }

    /// Few colors and sparse edges keep ind[S] near 2^n, exercising the counter width up to n = 24
    for (long long it = 0; it < stress::scaled(6); it++){
        int n = stress::rand_int(17, 24), c = stress::rand_int(1, 3);
        assert(planted(n, c) == c);
    }

    assert(planted(24, 4) == 4);
    return 0;
}
