/***
 *
 * Kruskal Reconstruction Tree
 * The vertices reachable from v using only edges of weight <= w, as one contiguous range of a fixed vertex order
 *
 * Complexity: O(m log m + n log n) to build, O(log n) per query, O(n log n) memory
 *
 * Kruskal's algorithm creates a new node for every merge, holding the merging edge's weight with the two merged
 * components as children: vertices 0..n-1 are the leaves, merges get ids n, n + 1, ..., 2n - 2, so a connected graph has 2n - 1 nodes
 * Weights never decrease going up, so v's component under threshold w is the subtree of v's highest ancestor with
 * weight <= w (found by binary lifting), and the leaves of every subtree are contiguous in order
 *
 * Usage:
 *   KruskalReconstructionTree<long long> krt(n); krt.add_edge(u, v, w); ... krt.build();   0-based vertices
 *   krt.range(v, w): half-open [l, r) of positions in krt.order holding exactly v's component
 *   krt.size(v, w): r - l, the number of vertices in v's component
 *   krt.node(v, w): the tree node of the component, v itself when it is alone
 *   krt.weight[x] for a merge node x >= n: the weight of the edge that created it, the component's bottleneck
 *   krt.order[i]: the vertex at position i, krt.pos[v]: the position of vertex v
 *   Add every edge before build, disconnected graphs, self-loops, multi-edges and negative weights are fine
 *   A segment tree over krt.order then aggregates "everything reachable from v with edges <= w"
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T = long long>
struct KruskalReconstructionTree{
    int n, lg = 1;
    vector<tuple<T, int, int>> edges;
    vector<T> weight;
    vector<int> order, pos, first, leaf_count;
    vector<vector<int>> up;

    KruskalReconstructionTree(int n) : n(n) {}

    void add_edge(int u, int v, T w){
        edges.push_back({w, u, v});
    }

    void build(){
        sort(edges.begin(), edges.end());
        vector<int> dsu(n), parent(n);
        vector<array<int, 2>> child(n, {-1, -1});
        iota(dsu.begin(), dsu.end(), 0);
        iota(parent.begin(), parent.end(), 0);
        weight.assign(n, T());

        auto find = [&](int x){
            int root = x;
            while (dsu[root] != root) root = dsu[root];
            while (dsu[x] != root){
                int next = dsu[x];
                dsu[x] = root, x = next;
            }
            return root;
        };

        for (auto [w, u, v] : edges){
            u = find(u), v = find(v);
            if (u == v) continue;
            int id = weight.size();
            dsu[u] = dsu[v] = parent[u] = parent[v] = id;
            dsu.push_back(id), parent.push_back(id), child.push_back({u, v}), weight.push_back(w);
        }
        int m = weight.size();

        /// children always have smaller ids than their parent, so id order is bottom-up and reverse id order top-down
        leaf_count.assign(m, 1);
        for (int x = n; x < m; x++) leaf_count[x] = leaf_count[child[x][0]] + leaf_count[child[x][1]];

        first.assign(m, -1);
        int offset = 0;
        for (int x = m - 1; x >= 0; x--){
            if (parent[x] == x) first[x] = offset, offset += leaf_count[x];
            if (x >= n) first[child[x][0]] = first[x], first[child[x][1]] = first[x] + leaf_count[child[x][0]];
        }

        order.assign(n, 0), pos.assign(n, 0);
        for (int v = 0; v < n; v++) pos[v] = first[v], order[first[v]] = v;

        while ((1 << lg) < m) lg++;
        up.assign(lg, parent);
        for (int k = 1; k < lg; k++){
            for (int x = 0; x < m; x++) up[k][x] = up[k - 1][up[k - 1][x]];
        }
    }

    /// a root is its own parent, so jumping past the root lands on the root, which is correct when its weight is <= w
    int node(int v, T w) const{
        for (int k = lg - 1; k >= 0; k--){
            if (weight[up[k][v]] <= w) v = up[k][v];
        }
        return v;
    }

    pair<int, int> range(int v, T w) const{
        int x = node(v, w);
        return {first[x], first[x] + leaf_count[x]};
    }

    int size(int v, T w) const{
        return leaf_count[node(v, w)];
    }
};

int main(){
    /***
     *   0 ==4,3== 1         3 --5-- 4 --1-- 5 (self-loop 0)
     *    \       /          |
     *     7     2           9
     *      \   /            |
     *        2 -------------+         6 isolated
    ***/
    KruskalReconstructionTree<int> krt(7);
    vector<array<int, 3>> edges = {{0, 1, 4}, {1, 2, 2}, {2, 0, 7}, {3, 4, 5}, {4, 5, 1}, {2, 3, 9}, {5, 5, 0}, {0, 1, 3}};
    for (auto [u, v, w] : edges) krt.add_edge(u, v, w);
    krt.build();

    auto component = [&](int v, int w){
        auto [l, r] = krt.range(v, w);
        vector<int> res(krt.order.begin() + l, krt.order.begin() + r);
        sort(res.begin(), res.end());
        assert(krt.size(v, w) == r - l);
        return res;
    };

    for (int v = 0; v < 7; v++) assert((component(v, 0) == vector<int>{v}));
    assert((component(5, 1) == vector<int>{4, 5}));
    assert((component(3, 1) == vector<int>{3}));
    assert((component(2, 2) == vector<int>{1, 2}));
    assert((component(0, 2) == vector<int>{0}));
    assert((component(0, 3) == vector<int>{0, 1, 2}));
    assert((component(3, 4) == vector<int>{3}));
    assert((component(3, 5) == vector<int>{3, 4, 5}));
    assert((component(1, 8) == vector<int>{0, 1, 2}));
    assert((component(4, 9) == vector<int>{0, 1, 2, 3, 4, 5}));
    assert((component(6, 1000) == vector<int>{6}));
    assert(krt.size(0, -5) == 1 && krt.size(2, 1000) == 6);
    assert(krt.weight[krt.node(4, 6)] == 5 && krt.weight[krt.node(0, 100)] == 9);

    for (int v = 0; v < 7; v++) assert(krt.order[krt.pos[v]] == v);

    KruskalReconstructionTree<long long> negative(3);
    negative.add_edge(0, 1, -1000000000000LL);
    negative.add_edge(1, 2, -5);
    negative.build();
    assert(negative.size(0, -1000000000001LL) == 1);
    assert(negative.size(0, -1000000000000LL) == 2);
    assert(negative.size(2, -6) == 1 && negative.size(2, -5) == 3);

    KruskalReconstructionTree<int> single(1);
    single.build();
    assert(single.size(0, 0) == 1 && single.range(0, 7) == make_pair(0, 1));

    return 0;
}
