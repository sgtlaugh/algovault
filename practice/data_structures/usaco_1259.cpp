/***
 *
 * Disjoint Set Union
 * Union-find: merges sets of nodes and tells whether two nodes share a set
 *
 * Complexity: amortized O(α(n)) per operation, union by size with path halving
 *
 * DSU dsu(n): nodes 0..n, each in its own set
 * dsu.connect(a, b): merge the sets of a and b, false if they are already one set
 * dsu.same(a, b): true if a and b are in the same set
 * dsu.size(a): number of nodes in a's set
 * dsu.find(a): the root of a's set, the same node for every member
 *
 * Variants: rollback_dsu.cpp (undoable unions), potential_dsu.cpp (relations x[b] = x[a] + d), partially_persistent_dsu.cpp (queries at past times)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DSU{
    vector<int> sz, parent;

    DSU(int n) : sz(n + 1, 1), parent(n + 1){
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int i){
        while (i != parent[i]){
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return parent[i];
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (sz[a] > sz[b]) swap(a, b);
        parent[a] = b, sz[b] += sz[a];
        return true;
    }

    bool same(int a, int b){
        return find(a) == find(b);
    }

    int size(int i){
        return sz[find(i)];
    }
};

int main(){
    int n, m;
    scanf("%d %d", &n, &m);
    vector<vector<int>> adj(n + 1);
    for (int i = 0; i < m; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u].push_back(v), adj[v].push_back(u);
    }

    vector<int> deg(n + 1), removed_deg(n + 1), order;
    vector<char> removed(n + 1, 0);
    set<pair<int, int>> alive;
    for (int i = 1; i <= n; i++) deg[i] = adj[i].size(), alive.insert({deg[i], i});
    while (!alive.empty()){
        int v = alive.begin()->second;
        alive.erase(alive.begin());
        removed[v] = 1, removed_deg[v] = deg[v], order.push_back(v);
        for (int u : adj[v]){
            if (removed[u]) continue;
            alive.erase({deg[u], u});
            alive.insert({--deg[u], u});
        }
    }

    /// v had the minimum degree when peeled, so every cow peeled after it keeps at least that many friends
    /// among those cows; the best group's first peeled cow sees the whole group in its component
    DSU dsu(n);
    vector<char> added(n + 1, 0);
    long long res = 0;
    for (int i = n - 1; i >= 0; i--){
        int v = order[i];
        added[v] = 1;
        for (int u : adj[v]){
            if (added[u]) dsu.connect(u, v);
        }
        res = max(res, (long long)removed_deg[v] * dsu.size(v));
    }
    printf("%lld\n", res);
    return 0;
}
