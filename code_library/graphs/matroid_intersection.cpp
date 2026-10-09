/***
 *
 * Matroid Intersection (unweighted)
 * Largest set of ground elements that is independent in two matroids at once
 *
 * Complexity: O(n * r^2 + n) oracle queries plus r + 1 oracle builds, n = ground set size, r = answer size
 *     with the oracles below every query is O(1), so O(n * r^2 + n + (r + 1) * build) overall
 *
 * matroid_intersection(n, m1, m2): ground elements 0..n - 1, returns the chosen elements sorted ascending
 * The oracles are template parameters, any struct with these three methods works:
 *     void build(const vector<int>& chosen)     chosen is the current common independent set, called before each round
 *     bool can_add(int x)                       is chosen + x independent, x not in chosen
 *     bool can_exchange(int y, int x)           is chosen - y + x independent, y in chosen, x not in chosen
 *
 * GraphicMatroid g(v): g.add_edge(a, b) adds element g.edges.size(), independent sets are forests
 *     on vertices 0..v - 1, self-loops are never independent, build is O(v + n)
 * PartitionMatroid p(color, cap): element x has color[x], at most cap[c] elements of color c, build is O(colors + r)
 *     colorful matroid: every cap is 1
 * XorMatroid x(vec): independent sets are linearly independent over GF(2), vectors are 64-bit, 0 is never
 *     independent, build is O(64 * n)
 *
 * Each round augments along a shortest path of the exchange graph, which keeps chosen independent in both
 * Bipartite matching is the intersection of two partition matroids, a rainbow spanning tree is graphic x colorful
 *
 * Example: largest forest with at most one edge per color
 *     GraphicMatroid g(v);  for each edge: g.add_edge(a, b)
 *     PartitionMatroid p(color, vector<int>(colors, 1));
 *     vector<int> chosen = matroid_intersection(m, g, p);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename M1, typename M2>
vector<int> matroid_intersection(int n, M1& m1, M2& m2){
    const int SOURCE = -1, UNSEEN = -2;
    vector<char> in_set(n, 0);
    vector<int> chosen;

    while (true){
        m1.build(chosen), m2.build(chosen);

        vector<int> parent(n, UNSEEN), bfs;
        for (int x = 0; x < n; x++){
            if (!in_set[x] && m1.can_add(x)) parent[x] = SOURCE, bfs.push_back(x);
        }

        int sink = -1;
        for (int i = 0; i < (int)bfs.size() && sink == -1; i++){
            int u = bfs[i];
            if (in_set[u]){
                for (int x = 0; x < n; x++){
                    if (!in_set[x] && parent[x] == UNSEEN && m1.can_exchange(u, x)) parent[x] = u, bfs.push_back(x);
                }
            }
            else if (m2.can_add(u)) sink = u;
            else{
                for (int y : chosen){
                    if (parent[y] == UNSEEN && m2.can_exchange(y, u)) parent[y] = u, bfs.push_back(y);
                }
            }
        }
        if (sink == -1) break;

        for (int x = sink; x != SOURCE; x = parent[x]) in_set[x] ^= 1;
        chosen.clear();
        for (int x = 0; x < n; x++){
            if (in_set[x]) chosen.push_back(x);
        }
    }
    return chosen;
}

/// chosen - y + x is a forest iff x joins two trees of chosen, or y lies on the tree path between x's endpoints
struct GraphicMatroid{
    int v;
    vector<array<int, 2>> edges;
    vector<vector<int>> adj;
    vector<int> comp, tin, tout, child, it;

    GraphicMatroid(int v) : v(v) {}

    void add_edge(int a, int b){
        edges.push_back({a, b});
    }

    void build(const vector<int>& chosen){
        adj.assign(v, {});
        for (int e : chosen){
            adj[edges[e][0]].push_back(e);
            adj[edges[e][1]].push_back(e);
        }

        comp.assign(v, -1), tin.assign(v, 0), tout.assign(v, 0), it.assign(v, 0);
        child.assign(edges.size(), -1);
        int timer = 0;
        for (int root = 0; root < v; root++){
            if (comp[root] != -1) continue;
            comp[root] = root, tin[root] = timer++;
            vector<int> stack = {root};
            while (!stack.empty()){
                int a = stack.back();
                if (it[a] == (int)adj[a].size()){
                    tout[a] = timer;
                    stack.pop_back();
                    continue;
                }

                int e = adj[a][it[a]++], b = edges[e][0] ^ edges[e][1] ^ a;
                if (comp[b] != -1) continue;
                comp[b] = root, tin[b] = timer++, child[e] = b;
                stack.push_back(b);
            }
        }
    }

    bool can_add(int x) const{
        return comp[edges[x][0]] != comp[edges[x][1]];
    }

    bool can_exchange(int y, int x) const{
        if (can_add(x)) return true;
        int c = child[y];
        return inside(c, edges[x][0]) != inside(c, edges[x][1]);
    }

    bool inside(int c, int a) const{
        return tin[c] <= tin[a] && tin[a] < tout[c];
    }
};

struct PartitionMatroid{
    vector<int> color, cap, used;

    PartitionMatroid(const vector<int>& color, const vector<int>& cap) : color(color), cap(cap) {}

    void build(const vector<int>& chosen){
        used.assign(cap.size(), 0);
        for (int x : chosen) used[color[x]]++;
    }

    bool can_add(int x) const{
        return used[color[x]] < cap[color[x]];
    }

    bool can_exchange(int y, int x) const{
        return color[x] == color[y] || can_add(x);
    }
};

/// Writes every element over the basis of chosen, then chosen - y + x is independent iff x is outside span(chosen) or uses y
struct XorMatroid{
    vector<unsigned long long> vec, uses;
    vector<char> spanned;
    vector<int> pos;

    XorMatroid(const vector<unsigned long long>& vec) : vec(vec) {}

    void build(const vector<int>& chosen){
        int n = vec.size();
        unsigned long long basis[64] = {}, basis_uses[64] = {};
        pos.assign(n, -1);
        for (int i = 0; i < (int)chosen.size(); i++){
            unsigned long long a = vec[chosen[i]], mask = 1ULL << i;
            pos[chosen[i]] = i;
            for (int b = 63; b >= 0; b--){
                if (!(a >> b & 1)) continue;
                if (!basis[b]){
                    basis[b] = a, basis_uses[b] = mask;
                    break;
                }
                a ^= basis[b], mask ^= basis_uses[b];
            }
        }

        uses.assign(n, 0), spanned.assign(n, 0);
        for (int x = 0; x < n; x++){
            unsigned long long a = vec[x], mask = 0;
            for (int b = 63; b >= 0; b--){
                if ((a >> b & 1) && basis[b]) a ^= basis[b], mask ^= basis_uses[b];
            }
            spanned[x] = a == 0, uses[x] = mask;
        }
    }

    bool can_add(int x) const{
        return !spanned[x];
    }

    bool can_exchange(int y, int x) const{
        return !spanned[x] || (uses[x] >> pos[y] & 1);
    }
};

int main(){
    GraphicMatroid tree(4);
    tree.add_edge(0, 1), tree.add_edge(1, 2), tree.add_edge(2, 3), tree.add_edge(0, 2), tree.add_edge(1, 3);
    PartitionMatroid one_per_color({0, 0, 1, 1, 2}, {1, 1, 1});
    vector<int> rainbow = matroid_intersection(5, tree, one_per_color);
    assert((rainbow == vector<int>{0, 2, 4} || rainbow == vector<int>{0, 3, 4} || rainbow == vector<int>{1, 3, 4}));

    GraphicMatroid parallel(4);
    parallel.add_edge(0, 1), parallel.add_edge(0, 1), parallel.add_edge(1, 0), parallel.add_edge(2, 3), parallel.add_edge(1, 2);
    PartitionMatroid colors({0, 1, 2, 0, 0}, {1, 1, 1});
    assert(matroid_intersection(5, parallel, colors).size() == 2);
    assert(matroid_intersection(5, colors, parallel).size() == 2);

    PartitionMatroid left({0, 0, 1, 2}, {1, 1, 1}), right({0, 1, 0, 0}, {1, 1});
    vector<int> matching = matroid_intersection(4, left, right);
    assert((matching == vector<int>{1, 2} || matching == vector<int>{1, 3}));

    XorMatroid vectors({1, 1, 2});
    PartitionMatroid swap_colors({1, 0, 1}, {1, 1});
    assert((matroid_intersection(3, swap_colors, vectors) == vector<int>{1, 2}));
    assert((matroid_intersection(3, vectors, swap_colors) == vector<int>{1, 2}));

    XorMatroid high_bits({1ULL << 63, (1ULL << 63) | 1, 1, 0});
    PartitionMatroid roomy({0, 0, 0, 0}, {4});
    assert(matroid_intersection(4, high_bits, roomy).size() == 2);

    GraphicMatroid loops(2);
    loops.add_edge(0, 0), loops.add_edge(1, 1);
    XorMatroid zeros({0, 0});
    assert(matroid_intersection(2, loops, zeros).empty());

    GraphicMatroid nothing(3);
    PartitionMatroid no_colors({}, {});
    assert(matroid_intersection(0, nothing, no_colors).empty());

    return 0;
}
