/***
 *
 * Partially Persistent DSU
 * Union-find that remembers its past: were a and b in the same set after the first t connect calls?
 *
 * Complexity: O(log n) per operation
 *
 * PartiallyPersistentDSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they already were one set
 * dsu.find(a, t): representative of a's set after the first t connect calls
 * dsu.same(a, b, t): true if a and b were in the same set after the first t connect calls
 * dsu.size(a, t): number of nodes in a's set after the first t connect calls
 * dsu.now: number of connect calls so far, failed ones included
 *
 * Notes:
 *   - t = 0 is the state before any call, any t >= now is the current state
 *   - The first call that joined a and b is the smallest t with same(a, b, t), found by binary search
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PartiallyPersistentDSU{
    int now = 0;
    vector<int> parent, merged_at;                /// merged_at[i]: call that gave i a parent, INT_MAX while i is a root
    vector<vector<pair<int, int>>> size_history;  /// (call, size) of root i, calls increasing

    PartiallyPersistentDSU(int n) : parent(n + 1), merged_at(n + 1, INT_MAX), size_history(n + 1, {{0, 1}}){
        iota(parent.begin(), parent.end(), 0);
    }

    bool connect(int a, int b){
        now++;
        a = find(a, now), b = find(b, now);
        if (a == b) return false;

        int size_a = size_history[a].back().second, size_b = size_history[b].back().second;
        if (size_a > size_b) swap(a, b);
        parent[a] = b, merged_at[a] = now;
        size_history[b].push_back({now, size_a + size_b});
        return true;
    }

    int find(int a, int t){
        assert(t >= 0);
        while (a != parent[a] && merged_at[a] <= t) a = parent[a];
        return a;
    }

    bool same(int a, int b, int t){
        return find(a, t) == find(b, t);
    }

    int size(int a, int t){
        auto& history = size_history[find(a, t)];
        return prev(upper_bound(history.begin(), history.end(), make_pair(t, INT_MAX)))->second;
    }
};

int main(){
    /// Four calls, the third one fails but still counts
    PartiallyPersistentDSU dsu(4);
    assert(dsu.connect(1, 2));
    assert(dsu.connect(3, 4));
    assert(!dsu.connect(2, 1));
    assert(dsu.connect(2, 4));
    assert(dsu.now == 4);
    assert(!dsu.same(1, 4, 3));
    assert(dsu.same(1, 4, 4));
    assert(dsu.size(4, 2) == 2);
    assert(dsu.size(4, 100) == 4);

    return 0;
}
