/***
 *
 * Prufer Code
 * Bijection between labeled trees on n >= 2 nodes (0-based) and sequences of n - 2 values in [0, n)
 * The code lists, for each removal of the smallest-labeled leaf, the neighbor of that leaf
 *
 * Complexity: O(n) for both encode and decode, O(n) memory
 *
 * Encode:
 *   PruferCode tree(n); tree.add_edge(u, v); vector<int> code = tree.encode();
 *   add_edge asserts both endpoints are in [0, n), encode asserts n >= 2 and that the n - 1 edges form a tree
 *   No adjacency lists or recursion: a leaf's only neighbor is the xor of its remaining neighbors
 *
 * Decode:
 *   vector<pair<int, int>> edges = prufer_decode(code);   n = code.size() + 2
 *   Returns the n - 1 edges as {child, parent} with the tree rooted at n - 1, asserts every value is in [0, n)
 *
 * Every code of length n - 2 decodes to a distinct tree, so there are n^(n - 2) labeled trees (Cayley)
 * The number of times v appears in the code is deg(v) - 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PruferCode{
    int n, added_edges = 0;
    vector<int> degree, neighbor_xor;

    PruferCode(int n) : n(n), degree(n), neighbor_xor(n) {}

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n);
        degree[u]++, degree[v]++;
        neighbor_xor[u] ^= v, neighbor_xor[v] ^= u;
        added_edges++;
    }

    /// Leaves are removed in label order, ptr only moves forward: a node that turns into a leaf
    /// below ptr is the smallest leaf right away, any other new leaf is found later by the scan
    /// Node n - 1 is never removed, and a removed leaf is never revisited since it sits at or below ptr
    vector<int> encode() const{
        assert(n >= 2 && added_edges == n - 1);
        vector<int> deg = degree, rest = neighbor_xor, code(n - 2);
        int ptr = 0;
        while (ptr < n && deg[ptr] != 1) ptr++;
        assert(ptr < n);

        int leaf = ptr;
        for (int i = 0; i < n - 2; i++){
            int next = rest[leaf];
            code[i] = next;
            rest[next] ^= leaf;
            if (--deg[next] == 1 && next < ptr){
                leaf = next;
                continue;
            }

            /// A graph with n - 1 edges that is not a tree has a cycle, which runs the scan out of leaves
            ptr++;
            while (ptr < n && deg[ptr] != 1) ptr++;
            assert(ptr < n);
            leaf = ptr;
        }
        return code;
    }
};

/// Same leaf order as encode, degree[v] - 1 counts the occurrences of v still ahead in the code
/// After n - 2 removals only the current leaf and n - 1 remain
vector<pair<int, int>> prufer_decode(const vector<int>& code){
    int n = code.size() + 2;
    vector<int> degree(n, 1);
    for (int v : code){
        assert(0 <= v && v < n);
        degree[v]++;
    }

    int ptr = 0;
    while (degree[ptr] != 1) ptr++;
    int leaf = ptr;
    vector<pair<int, int>> edges;
    edges.reserve(n - 1);
    for (int v : code){
        edges.push_back({leaf, v});
        if (--degree[v] == 1 && v < ptr){
            leaf = v;
            continue;
        }
        while (degree[++ptr] != 1);
        leaf = ptr;
    }

    edges.push_back({leaf, n - 1});
    return edges;
}

int main(){
    PruferCode star(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 3}, {1, 3}, {2, 3}, {3, 4}, {4, 5}}) star.add_edge(u, v);
    assert((star.encode() == vector<int>{3, 3, 3, 4}));
    assert((prufer_decode({3, 3, 3, 4}) == vector<pair<int, int>>{{0, 3}, {1, 3}, {2, 3}, {3, 4}, {4, 5}}));

    PruferCode path(4);
    path.add_edge(2, 3), path.add_edge(1, 2), path.add_edge(0, 1);
    assert((path.encode() == vector<int>{1, 2}));
    assert((prufer_decode({1, 2}) == vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 3}}));

    PruferCode early_leaf(5);
    early_leaf.add_edge(4, 3), early_leaf.add_edge(1, 0), early_leaf.add_edge(4, 2), early_leaf.add_edge(0, 4);
    assert((early_leaf.encode() == vector<int>{0, 4, 4}));
    assert((prufer_decode({0, 4, 4}) == vector<pair<int, int>>{{1, 0}, {0, 4}, {2, 4}, {3, 4}}));

    assert((prufer_decode({2, 2}) == vector<pair<int, int>>{{0, 2}, {1, 2}, {2, 3}}));
    assert((prufer_decode({0, 0, 0}) == vector<pair<int, int>>{{1, 0}, {2, 0}, {3, 0}, {0, 4}}));

    PruferCode pair_tree(2);
    pair_tree.add_edge(1, 0);
    assert(pair_tree.encode().empty());
    assert((prufer_decode({}) == vector<pair<int, int>>{{0, 1}}));

    return 0;
}
