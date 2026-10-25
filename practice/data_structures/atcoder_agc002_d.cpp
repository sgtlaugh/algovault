/***
 *
 * Partially Persistent DSU
 * Union-find that remembers its past: were a and b in the same set after the first t connect calls?
 *
 * Complexity: O(log n) per operation
 *
 * PartiallyPersistentDSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they already were one set
 * dsu.same(a, b, t): true if a and b were in the same set after the first t connect calls
 * dsu.size(a, t): number of nodes in a's set after the first t connect calls
 * dsu.find(a, t): representative of a's set after the first t connect calls
 * dsu.now: number of connect calls so far, failed ones included
 *
 * t = 0 is the state before any call, any t >= now is the current state
 * The first call that joined a and b is the smallest t with same(a, b, t), found by binary search
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

    int find(int a, int t){
        assert(t >= 0);
        while (a != parent[a] && merged_at[a] <= t) a = parent[a];
        return a;
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

    bool same(int a, int b, int t){
        return find(a, t) == find(b, t);
    }

    int size(int a, int t){
        auto& history = size_history[find(a, t)];
        /// the last entry recorded at or before call t
        return prev(upper_bound(history.begin(), history.end(), make_pair(t, INT_MAX)))->second;
    }
};

int main(){
    int n, m;

    scanf("%d %d", &n, &m);
    PartiallyPersistentDSU dsu(n);
    for (int i = 0; i < m; i++){
        int a, b;
        scanf("%d %d", &a, &b);
        dsu.connect(a, b);
    }

    auto visited = [&](int x, int y, int t){
        return dsu.same(x, y, t) ? dsu.size(x, t) : dsu.size(x, t) + dsu.size(y, t);
    };

    int q;
    scanf("%d", &q);
    while (q--){
        int x, y, z;
        scanf("%d %d %d", &x, &y, &z);
        int lo = 1, hi = m;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (visited(x, y, mid) >= z) hi = mid;
            else lo = mid + 1;
        }
        printf("%d\n", lo);
    }
    return 0;
}
