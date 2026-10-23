/***
 *
 * Prufer Code
 * Bijection between labeled trees on n >= 2 nodes (0-based) and sequences of n - 2 values in [0, n)
 *
 * Complexity: O(n) for both encode and decode, O(n) memory
 *
 * PruferCode tree(n): an empty tree on nodes 0..n-1
 * tree.add_edge(u, v): adds the edge u - v, asserts both endpoints are in [0, n)
 * tree.encode(): for each removal of the smallest-labeled leaf, its neighbor, asserts n >= 2 and a tree
 * prufer_decode(code): the tree on n = code.size() + 2 nodes as n - 1 edges {child, parent}, rooted at n - 1
 *
 * Notes:
 *   - prufer_decode asserts every value is in [0, n)
 *   - v appears deg(v) - 1 times in the code
 *   - Every code decodes to a distinct tree, so there are n^(n - 2) labeled trees (Cayley)
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

    /// Leaves go in label order with ptr only moving forward: a new leaf below ptr is the smallest leaf right away
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
    /***
     *   0   1   2
     *    \  |  /
     *       3 - 4 - 5
    ***/
    PruferCode tree(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 3}, {1, 3}, {2, 3}, {3, 4}, {4, 5}}) tree.add_edge(u, v);
    assert((tree.encode() == vector<int>{3, 3, 3, 4}));  /// leaves 0, 1, 2 hang off 3, then 3 itself is a leaf of 4

    vector<pair<int, int>> edges = prufer_decode({3, 3, 3, 4});
    assert((edges == vector<pair<int, int>>{{0, 3}, {1, 3}, {2, 3}, {3, 4}, {4, 5}}));  /// {child, parent}, rooted at 5
    return 0;
}
