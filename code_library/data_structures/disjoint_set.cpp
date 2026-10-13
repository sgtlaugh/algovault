/***
 *
 * Disjoint Set Union
 * Union by size with path compression
 * Nodes can be 0 or 1 based, allocates nodes 0..n
 * connect returns true if a and b were in different components
 *
 * Complexity: amortized O(α(n)) per operation, where α(n) is the inverse Ackermann function, O(n) memory
 *
 * Variants: rollback_dsu.cpp (undoable unions), potential_dsu.cpp (relations x[b] = x[a] + d), partially_persistent_dsu.cpp (queries at past times)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DSU{
    vector<int> sz, parent;

    DSU(int n){
        parent.resize(n + 1);
        sz.resize(n + 1, 1);
        for (int i = 0; i <= n; i++) parent[i] = i;
    }

    int find(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return parent[i];
    }

    bool connect(int a, int b){
        int c = find(a), d = find(b);
        if (c == d) return false;
        if (sz[c] > sz[d]) swap(c, d);
        parent[c] = d;
        sz[d] += sz[c], sz[c] = 0;
        return true;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        int p = find(i);
        return sz[p];
    }
};

int main(){
    auto dsu = DSU(8);
    dsu.connect(1, 2);
    dsu.connect(1, 3);
    dsu.connect(4, 5);
    dsu.connect(5, 6);
    dsu.connect(6, 7);

    assert(dsu.size(3) == 3);
    assert(dsu.size(5) == 4);
    assert(!dsu.same(3, 5));

    assert(dsu.connect(2, 5));
    assert(!dsu.connect(7, 1));

    assert(dsu.size(2) == 7);
    assert(dsu.size(6) == 7);
    assert(dsu.same(3, 5));

    return 0;
}
