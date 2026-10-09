/***
 *
 * Edge Coloring
 * Colors the edges so that edges sharing a vertex get different colors
 *
 * Complexity: BipartiteEdgeColoring O(n + m * n) expected time, O(n + m) memory
 * VizingEdgeColoring O(n * m) time, O(n * D) memory
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

        /// A free color at x always exists below deg(x): candidates[s][x] holds every free one there, stale ones are popped lazily
        vector<vector<int>> candidates[2];
        for (int s = 0; s < 2; s++){
            for (int d : deg[s]){
                candidates[s].emplace_back(d);
                iota(candidates[s].back().rbegin(), candidates[s].back().rend(), 0);
            }
        }

        /// Hash map of (side, vertex, color) to the edge colored so there, -1 if none, instead of an O(n * D) table
        /// At most 4m keys ever appear (colors below the degree, plus two per edge), so 8m slots keep probing short
        int bits = 1;
        while ((1LL << bits) < 8LL * m) bits++;
        vector<pair<long long, int>> table(1LL << bits, {-1, -1});
        vector<int> color(m, -1), path;
        auto slot = [&](int s, int x, int c) -> int& {
            long long key = (s ? (long long)n_left + x : x) * k + c;
            size_t i = (unsigned long long)key * 0x9E3779B97F4A7C15ULL >> (64 - bits);
            while (table[i].first != key && table[i].first != -1) i = (i + 1) & (table.size() - 1);
            if (table[i].first == -1) table[i] = {key, -1};
            return table[i].second;
        };
        auto first_free = [&](int s, int x){
            auto& list = candidates[s][x];
            while (slot(s, x, list.back()) != -1) list.pop_back();
            return list.back();
        };

        for (int e = 0; e < m; e++){
            int u = edges[e][0], v = edges[e][1], a = first_free(0, u), b = first_free(1, v);
            if (slot(1, v, a) != -1){
                /// The a/b path from v enters left vertices by a edges, so it never reaches u where a is free
                path.clear();
                int side = 1, x = v, c = a;
                for (; slot(side, x, c) != -1; side ^= 1, c ^= a ^ b){
                    int f = slot(side, x, c);
                    path.push_back(f);
                    x = edges[f][side ^ 1];
                }
                /// The path's far end loses color c ^ a ^ b, the only color the inversion frees apart from a at v
                candidates[side][x].push_back(c ^ a ^ b);

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
        vector<int> adj((size_t)n * k, -1), free_color(n), fan(k + 1), fan_color(k + 1), seen;
        auto at = [&](int x, int c) -> int& { return adj[(long long)x * k + c]; };

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
    BipartiteEdgeColoring path(2, 2);
    path.add_edge(0, 0), path.add_edge(1, 0), path.add_edge(1, 1);
    vector<int> color = path.solve();
    assert(color[0] != color[1]);   /// both touch right vertex 0
    assert(color[1] != color[2]);   /// both touch left vertex 1
    assert(color[0] == color[2]);   /// max degree 2, so only colors 0 and 1 exist

    /// A triangle needs D + 1 = 3 colors, the most Vizing ever uses
    VizingEdgeColoring triangle(3);
    triangle.add_edge(0, 1), triangle.add_edge(1, 2), triangle.add_edge(2, 0);
    color = triangle.solve();
    sort(color.begin(), color.end());
    assert((color == vector<int>{0, 1, 2}));
    return 0;
}
