/***
 *
 * Leaf Peeling
 * Rooted parents and a leaves-first order of a tree given as an edge list, with no adjacency lists
 *
 * Complexity: O(n), two ints per node
 *
 * PeeledTree t = peel_leaves(n, edges, root = 0): edges are the n - 1 pairs {u, v} of a tree on nodes 0 to n - 1
 * t.order: all n nodes, every node after all of its children, root last
 * t.parent: parent[v] when the tree hangs from root, -1 for root
 *
 * Notes:
 *   - Bottom-up DP such as subtree sizes walks order; top-down DP such as depths walks it backwards
 *   - Children are never listed and order is not a DFS order, so use an adjacency list when either is needed
 *   - Asserts that the edges form a tree
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PeeledTree{
    vector<int> order, parent;
};

/// rest[v] is the XOR of v's neighbors not yet peeled, so a leaf's rest is its parent
PeeledTree peel_leaves(int n, const vector<pair<int, int>>& edges, int root = 0){
    assert(n >= 1 && (int)edges.size() == n - 1 && 0 <= root && root < n);
    vector<int> deg(n, 0), rest(n, 0);
    for (auto [u, v] : edges){
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        deg[u]++, deg[v]++;
        rest[u] ^= v, rest[v] ^= u;
    }

    PeeledTree res{{}, vector<int>(n, -1)};
    res.order.reserve(n);
    for (int v = 0; v < n; v++){
        if (deg[v] == 1 && v != root) res.order.push_back(v);
    }

    for (int i = 0; i < (int)res.order.size(); i++){
        int v = res.order[i], p = rest[v];
        assert(deg[v] == 1);
        res.parent[v] = p;
        rest[p] ^= v;
        if (--deg[p] == 1 && p != root) res.order.push_back(p);
    }

    res.order.push_back(root);
    assert((int)res.order.size() == n);
    return res;
}

int main(){
    /// Rooted at 0: 1 and 2 hang from 0, 3 and 4 from 1; leaves 2, 3, 4 come first, then 1, then 0
    PeeledTree t = peel_leaves(5, {{0, 1}, {1, 3}, {0, 2}, {4, 1}});
    assert((t.order == vector<int>{2, 3, 4, 1, 0}));
    assert((t.parent == vector<int>{-1, 0, 0, 1, 1}));

    vector<int> size(5, 1);
    for (int v : t.order){
        if (t.parent[v] != -1) size[t.parent[v]] += size[v];
    }
    assert(size[1] == 3);

    PeeledTree from_three = peel_leaves(5, {{0, 1}, {1, 3}, {0, 2}, {4, 1}}, 3);
    assert((from_three.parent == vector<int>{1, 3, 0, -1, 1}));

    return 0;
}
