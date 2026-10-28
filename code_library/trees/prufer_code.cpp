/***
 *
 * Prufer Code
 * Bijection between labeled trees on n >= 2 nodes (0-based) and sequences of n - 2 values in [0, n)
 *
 * Complexity: O(n) for both encode and decode
 *
 * prufer_encode(n, edges): for each removal of the smallest-labeled leaf, its neighbor, asserts edges form a tree
 * prufer_decode(code): the tree on n = code.size() + 2 nodes as n - 1 edges {child, parent}, rooted at n - 1
 *
 * Notes:
 *   - v appears deg(v) - 1 times in the code
 *   - Every code decodes to a distinct tree, so there are n^(n - 2) labeled trees (Cayley)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Leaves go in label order with ptr only moving forward: a new leaf below ptr is the smallest leaf right away
vector<int> prufer_encode(int n, const vector<pair<int, int>>& edges){
    assert(n >= 2 && (int)edges.size() == n - 1);
    vector<int> deg(n), rest(n), code(n - 2);
    for (auto [u, v] : edges){
        assert(0 <= u && u < n && 0 <= v && v < n);
        deg[u]++, deg[v]++;
        rest[u] ^= v, rest[v] ^= u;
    }

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

/// Same leaf order as prufer_encode, degree[v] - 1 counts the occurrences of v still ahead in the code
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
    vector<pair<int, int>> edges = {{0, 3}, {1, 3}, {2, 3}, {3, 4}, {4, 5}};
    assert((prufer_encode(6, edges) == vector<int>{3, 3, 3, 4}));  /// leaves 0, 1, 2 hang off 3, then 3 is a leaf of 4

    assert(prufer_decode({3, 3, 3, 4}) == edges);  /// the same tree as {child, parent} edges, rooted at 5
    return 0;
}
