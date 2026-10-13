/***
 *
 * Disjoint Set Union With Rollback
 * Undoable unions, for offline dynamic connectivity and divide and conquer over time
 * Nodes can be 0 or 1 based, allocates nodes 0..n
 * connect returns true if a and b were in different components
 *
 * Complexity: O(log n) per operation, union by size without path compression, O(n) memory
 *
 * time(): number of successful connects not yet rolled back
 * rollback(t): undo every union made after time() returned t, O(1) per undone union
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
    auto rb = RollbackDSU(6);
    assert(rb.connect(1, 2));
    int checkpoint = rb.time();
    assert(rb.connect(3, 4) && rb.connect(2, 3) && !rb.connect(1, 4));

    assert(rb.time() == 3 && checkpoint == 1);
    assert(rb.size(4) == 4 && rb.same(1, 4));

    rb.rollback(checkpoint);
    assert(rb.time() == 1);
    assert(rb.size(1) == 2 && rb.size(3) == 1 && rb.size(4) == 1);
    assert(!rb.same(1, 3) && !rb.same(3, 4));

    assert(rb.connect(4, 5) && rb.size(5) == 2);
    rb.rollback(0);
    for (int i = 0; i <= 6; i++) assert(rb.size(i) == 1);

    return 0;
}
