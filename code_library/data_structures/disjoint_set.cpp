/***
 *
 * Disjoint Set Union
 * Union-find: merges sets of nodes and tells whether two nodes share a set
 *
 * Complexity: amortized O(α(n)) per operation
 *
 * DSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they are already one set
 * dsu.find(a): the root of a's set, the same node for every member
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 *
 * Notes:
 *   - Use rollback_dsu.cpp to undo unions, potential_dsu.cpp for relations x[b] = x[a] + d
 *   - Use partially_persistent_dsu.cpp for queries at past times
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DSU{
    vector<int> sz, parent;

    DSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        return true;
    }

    int find(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }

        return parent[i];
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }
};

int main(){
    /// Groups {1, 2, 3} and {4, 5, 6, 7} merge into one
    DSU dsu(8);
    assert(dsu.connect(1, 2));
    assert(dsu.connect(1, 3));
    assert(dsu.connect(4, 5));
    assert(dsu.connect(5, 6));
    assert(dsu.connect(6, 7));
    assert(dsu.size(3) == 3);
    assert(dsu.size(5) == 4);
    assert(!dsu.same(3, 5));

    assert(dsu.connect(2, 5));
    assert(!dsu.connect(7, 1));
    assert(dsu.size(6) == 7);
    assert(dsu.same(3, 5));
    assert(dsu.find(8) == 8);

    return 0;
}
