/***
 *
 * Edge Coloring
 * Colors the edges so that edges sharing a vertex get different colors
 *
 * Complexity: BipartiteEdgeColoring O(m * (n + D)), VizingEdgeColoring O(n * m), both O(n * D) memory
 * where n is the number of vertices, m the number of edges and D the maximum degree
 *
 * BipartiteEdgeColoring g(n_left, n_right): left vertices 0..n_left - 1, right vertices 0..n_right - 1
 * g.add_edge(u, v): edge from left u to right v, parallel edges allowed
 * g.solve(): color of every edge in insertion order, in 0..D - 1, optimal by Konig's theorem
 *
 * VizingEdgeColoring g(n): vertices 0..n - 1, the graph must be simple (no loops, no parallel edges)
 * g.add_edge(u, v): undirected edge
 * g.solve(): color of every edge in insertion order, in 0..D (Misra-Gries)
 * Finding a D coloring of a general graph is NP-hard, so D + 1 colors is the guarantee
 *
 * Both solve() calls can be repeated after adding more edges
 *
 * Example:
 *   BipartiteEdgeColoring g(2, 2);
 *   g.add_edge(0, 0), g.add_edge(0, 1), g.add_edge(1, 0);
 *   vector<int> color = g.solve();  // color[0] != color[1], color[0] != color[2], all in {0, 1}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct BipartiteEdgeColoring{
    int n_left, n_right;
    vector<array<int, 2>> edges;

    BipartiteEdgeColoring(int n_left, int n_right) : n_left(n_left), n_right(n_right) {}

    void add_edge(int u, int v){
        edges.push_back({u, v});
    }

    vector<int> solve(){
        int m = edges.size(), k = 0;
        vector<int> deg[2] = {vector<int>(n_left), vector<int>(n_right)};
        for (auto& e : edges){
            k = max({k, ++deg[0][e[0]], ++deg[1][e[1]]});
        }

        /// edge_at[s][x * k + c]: the edge colored c at vertex x of side s, or -1
        vector<int> edge_at[2] = {vector<int>(n_left * k, -1), vector<int>(n_right * k, -1)};
        vector<int> color(m, -1), path;
        auto slot = [&](int s, int x, int c) -> int& { return edge_at[s][x * k + c]; };
        auto first_free = [&](int s, int x){
            int c = 0;
            while (slot(s, x, c) != -1) c++;
            return c;
        };

        for (int e = 0; e < m; e++){
            int u = edges[e][0], v = edges[e][1], a = first_free(0, u), b = first_free(1, v);
            if (slot(1, v, a) != -1){
                /// The a/b path from v enters left vertices by a edges, so it never reaches u where a is free
                path.clear();
                for (int s = 1, x = v, c = a; slot(s, x, c) != -1; s ^= 1, c ^= a ^ b){
                    int f = slot(s, x, c);
                    path.push_back(f);
                    x = edges[f][s ^ 1];
                }
                for (int f : path){
                    for (int s = 0; s < 2; s++) slot(s, edges[f][s], color[f]) = -1;
                }
                for (int f : path){
                    color[f] ^= a ^ b;
                    for (int s = 0; s < 2; s++) slot(s, edges[f][s], color[f]) = f;
                }
            }
            color[e] = a;
            slot(0, u, a) = slot(1, v, a) = e;
        }
        return color;
    }
};

struct VizingEdgeColoring{
    int n;
    vector<array<int, 2>> edges;

    VizingEdgeColoring(int n) : n(n) {}

    void add_edge(int u, int v){
        assert(u != v);
        edges.push_back({u, v});
    }

    vector<int> solve(){
        int m = edges.size(), k = 1;
        vector<int> deg(n);
        for (auto& e : edges){
            k = max({k, ++deg[e[0]] + 1, ++deg[e[1]] + 1});
        }

        /// adj[x * k + c]: the neighbor of x through the edge colored c, or -1
        vector<int> adj(n * k, -1), free_color(n), fan(k + 1), fan_color(k + 1), seen;
        auto at = [&](int x, int c) -> int& { return adj[x * k + c]; };

        for (auto& e : edges){
            int u = e[0], v = e[1], c = free_color[u], d, len = 0;
            fan[0] = v;
            seen.assign(k, 0);
            while (d = free_color[v], !seen[d] && (v = at(u, d)) != -1){
                seen[d] = ++len;
                fan_color[len] = d, fan[len] = v;
            }

            /// Inverts the c/d path starting at u, which recolors u's d edge (fan[seen[d]], if any) to c
            fan_color[seen[d]] = c;
            int end = u;
            for (int x = u, cd = d; x != -1; cd ^= c ^ d, x = at(x, cd)){
                swap(at(x, cd), at(x, cd ^ c ^ d));
                end = x;
            }

            /// Rotates the fan until its current tip has d free, then colors the tip's edge d
            int i = 0;
            while (at(fan[i], d) != -1){
                int left = fan[i], right = fan[++i], f = fan_color[i];
                at(u, f) = left, at(left, f) = u;
                at(right, f) = -1;
                free_color[right] = f;
            }
            at(u, d) = fan[i], at(fan[i], d) = u;

            for (int y : {fan[0], u, end}){
                free_color[y] = 0;
                while (at(y, free_color[y]) != -1) free_color[y]++;
            }
        }

        vector<int> color(m, 0);
        for (int i = 0; i < m; i++){
            while (at(edges[i][0], color[i]) != edges[i][1]) color[i]++;
        }
        return color;
    }
};

int main(){
    auto is_proper = [](int n, const vector<array<int, 2>>& edges, const vector<int>& color, int num_colors){
        if (color.size() != edges.size()) return false;
        set<pair<int, int>> used;
        for (int i = 0; i < (int)edges.size(); i++){
            if (color[i] < 0 || color[i] >= num_colors) return false;
            for (int x : edges[i]){
                if (x < 0 || x >= n || !used.insert({x, color[i]}).second) return false;
            }
        }
        return true;
    };

    vector<array<int, 2>> k33;
    BipartiteEdgeColoring full(3, 3);
    for (int u = 0; u < 3; u++){
        for (int v = 0; v < 3; v++) full.add_edge(u, v), k33.push_back({u, 3 + v});
    }
    assert(is_proper(6, k33, full.solve(), 3));

    BipartiteEdgeColoring path(2, 2);
    path.add_edge(0, 0), path.add_edge(1, 0), path.add_edge(1, 1);
    vector<int> path_color = path.solve();
    assert(path_color[0] != path_color[1] && path_color[1] != path_color[2] && path_color[0] == path_color[2]);
    assert(path_color[0] < 2 && path_color[1] < 2);

    BipartiteEdgeColoring parallel(1, 1);
    for (int i = 0; i < 3; i++) parallel.add_edge(0, 0);
    vector<int> parallel_color = parallel.solve();
    sort(parallel_color.begin(), parallel_color.end());
    assert((parallel_color == vector<int>{0, 1, 2}));

    BipartiteEdgeColoring grow(2, 2);
    grow.add_edge(0, 0), grow.add_edge(1, 1);
    assert((grow.solve() == vector<int>{0, 0}));
    grow.add_edge(0, 1);
    assert(is_proper(4, {{0, 2}, {1, 3}, {0, 3}}, grow.solve(), 2));

    assert(BipartiteEdgeColoring(0, 0).solve().empty());
    assert(BipartiteEdgeColoring(4, 2).solve().empty());

    VizingEdgeColoring triangle(3);
    triangle.add_edge(0, 1), triangle.add_edge(1, 2), triangle.add_edge(2, 0);
    vector<int> triangle_color = triangle.solve();
    sort(triangle_color.begin(), triangle_color.end());
    assert((triangle_color == vector<int>{0, 1, 2}));

    vector<array<int, 2>> petersen;
    for (int i = 0; i < 5; i++){
        petersen.push_back({i, (i + 1) % 5});
        petersen.push_back({i, i + 5});
        petersen.push_back({i + 5, (i + 2) % 5 + 5});
    }
    VizingEdgeColoring snark(10);
    for (auto& e : petersen) snark.add_edge(e[0], e[1]);
    vector<int> petersen_color = snark.solve();
    assert(is_proper(10, petersen, petersen_color, 4));
    assert(*max_element(petersen_color.begin(), petersen_color.end()) == 3);

    vector<array<int, 2>> k4;
    VizingEdgeColoring clique(4);
    for (int u = 0; u < 4; u++){
        for (int v = u + 1; v < 4; v++) clique.add_edge(u, v), k4.push_back({u, v});
    }
    assert(is_proper(4, k4, clique.solve(), 4));

    VizingEdgeColoring star(5);
    for (int v = 1; v < 5; v++) star.add_edge(0, v);
    assert(is_proper(5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}}, star.solve(), 5));

    VizingEdgeColoring regrow(4);
    regrow.add_edge(0, 1);
    assert((regrow.solve() == vector<int>{0}));
    regrow.add_edge(1, 2), regrow.add_edge(2, 0), regrow.add_edge(0, 3);
    assert(is_proper(4, {{0, 1}, {1, 2}, {2, 0}, {0, 3}}, regrow.solve(), 4));

    VizingEdgeColoring single(2);
    single.add_edge(1, 0);
    assert(is_proper(2, {{1, 0}}, single.solve(), 2));

    assert(VizingEdgeColoring(0).solve().empty());
    assert(VizingEdgeColoring(6).solve().empty());

    return 0;
}
