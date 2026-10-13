/***
 *
 * Disjoint Set Union With Rollback
 * Union-find whose unions can be undone, for offline dynamic connectivity and divide and conquer over time
 *
 * Complexity: O(log n) per operation, O(1) per undone union
 *
 * RollbackDSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they are already one set
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 * dsu.find(a): the root of a's set, the same node for every member
 * dsu.time(): number of successful connects not yet rolled back, use it as a checkpoint
 * dsu.rollback(t): undo every union made after time() returned t, 0 <= t <= time()
 *
 * Union by size without path compression, so a union changes one parent and undoing it is O(1)
 * offline_dynamic_connectivity.cpp carries a checked copy of this struct
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED rollback_dsu
struct RollbackDSU{
    vector<int> sz, parent, history;

    RollbackDSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int i){
        while (i != parent[i]) i = parent[i];
        return i;
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        history.push_back(a);
        return true;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }

    int time(){
        return history.size();
    }

    void rollback(int t){
        assert(0 <= t && t <= time());
        while ((int)history.size() > t){
            int a = history.back();
            history.pop_back();
            sz[parent[a]] -= sz[a], parent[a] = a;
        }
    }
};
/// END SHARED rollback_dsu

int main(){
    /// Take a checkpoint, merge more, roll back to it
    RollbackDSU dsu(6);
    assert(dsu.connect(1, 2));
    int checkpoint = dsu.time();
    assert(checkpoint == 1);
    assert(dsu.connect(3, 4));
    assert(dsu.connect(2, 3));
    assert(!dsu.connect(1, 4));
    assert(dsu.time() == 3);
    assert(dsu.size(4) == 4);
    assert(dsu.same(1, 4));

    dsu.rollback(checkpoint);
    assert(dsu.time() == 1);
    assert(dsu.size(1) == 2);
    assert(dsu.size(4) == 1);
    assert(!dsu.same(1, 3));

    return 0;
}
