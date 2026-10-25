/***
 *
 * Weighted Matroid Intersection
 * Maximum weight set of ground elements that is independent in two matroids at once
 *
 * Complexity:
 *   - weighted_matroid_intersection: O(n r^2) oracle queries, r + 1 builds, O(n^2 r^2) time, for an answer of size r
 *   - GraphicMatroid: build O(v + n), queries O(1)
 *   - PartitionMatroid: build O(colors + r), queries O(1)
 *   - XorMatroid: build O(64 n), queries O(1)
 *
 * weighted_matroid_intersection(n, weight, m1, m2): a max weight set of elements 0 to n - 1 independent in both, sorted
 * m1 and m2 are oracles, any structs with these three calls:
 *     void build(const vector<int>& chosen): chosen is the current common independent set, called before each round
 *     bool can_add(int x): true if chosen + x is independent, x not in chosen
 *     bool can_exchange(int y, int x): true if chosen - y + x is independent, y in chosen, x not in chosen
 *
 * GraphicMatroid g(v): independent sets are forests on vertices 0 to v - 1, self-loops are never independent
 * g.add_edge(a, b): adds the edge a - b as element g.edges.size()
 * PartitionMatroid p(color, cap): element x has color[x], an independent set has at most cap[c] elements of color c
 * XorMatroid x(vec): independent sets of 64-bit vectors are linearly independent over GF(2), 0 never is
 *
 * Notes:
 *   - weight[x] may be negative, the sum of all |weight[x]| must fit in long long
 *   - Best among the largest sets: add B = 2 sum |weight[x]| + 1 to every weight, n (B + max |weight[x]|) must fit
 *   - Cheapest instead of heaviest: negate every weight
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Augments along the cheapest exchange path, ties broken by fewest elements, which keeps chosen independent in both
template<typename M1, typename M2>
vector<int> weighted_matroid_intersection(int n, const vector<long long>& weight, M1& m1, M2& m2){
    vector<char> in_set(n, 0);
    vector<int> chosen;

    while (true){
        m1.build(chosen), m2.build(chosen);

        vector<vector<int>> next(n);
        for (int y : chosen){
            for (int x = 0; x < n; x++){
                if (in_set[x]) continue;
                if (m1.can_exchange(y, x)) next[y].push_back(x);
                if (m2.can_exchange(y, x)) next[x].push_back(y);
            }
        }

        vector<long long> cost(n, LLONG_MAX);
        vector<int> hops(n, 0), parent(n, -1);
        for (int x = 0; x < n; x++){
            if (!in_set[x] && m1.can_add(x)) cost[x] = -weight[x];
        }
        for (bool changed = true; changed;){
            changed = false;
            for (int u = 0; u < n; u++){
                if (cost[u] == LLONG_MAX) continue;
                for (int v : next[u]){
                    long long c = cost[u] + (in_set[v] ? weight[v] : -weight[v]);
                    if (c < cost[v] || (c == cost[v] && hops[u] + 1 < hops[v])){
                        cost[v] = c, hops[v] = hops[u] + 1, parent[v] = u, changed = true;
                    }
                }
            }
        }

        int sink = -1;
        for (int x = 0; x < n; x++){
            if (in_set[x] || cost[x] == LLONG_MAX || !m2.can_add(x)) continue;
            if (sink == -1 || cost[x] < cost[sink] || (cost[x] == cost[sink] && hops[x] < hops[sink])) sink = x;
        }
        /// The best weight of a k element set is concave in k, so the first path that gains nothing ends the search
        if (sink == -1 || cost[sink] >= 0) break;

        for (int x = sink; x != -1; x = parent[x]) in_set[x] ^= 1;
        chosen.clear();
        for (int x = 0; x < n; x++){
            if (in_set[x]) chosen.push_back(x);
        }
    }

    return chosen;
}

/// BEGIN COPY matroid_oracles from code_library/graphs/matroid_intersection.cpp
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

/// Writes each element over the basis of chosen: chosen - y + x is independent iff x is outside span(chosen) or uses y
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
/// END COPY matroid_oracles

int main(){
    /// Edges 0-1 and 1-2 share color 0, so the heaviest rainbow forest is 0-1 and 0-2 at 5 + 3
    GraphicMatroid forest(3);
    forest.add_edge(0, 1);
    forest.add_edge(1, 2);
    forest.add_edge(0, 2);
    PartitionMatroid one_per_color({0, 0, 1}, {1, 1});
    assert((weighted_matroid_intersection(3, {5, 4, 3}, forest, one_per_color) == vector<int>{0, 2}));

    /// Max weight bipartite matching, element x joins left[x] to right[x]: the 3 blocks both 2s, which weigh 4 together
    PartitionMatroid left({0, 0, 1}, {1, 1});
    PartitionMatroid right({0, 1, 1}, {1, 1});
    assert((weighted_matroid_intersection(3, {2, 3, 2}, left, right) == vector<int>{0, 2}));

    return 0;
}
