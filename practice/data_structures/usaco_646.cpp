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
    freopen("closing.in", "r", stdin);
    freopen("closing.out", "w", stdout);

    int n, m;
    scanf("%d %d", &n, &m);

    vector<vector<int>> adj(n + 1);
    for (int i = 0; i < m; i++){
        int a, b;
        scanf("%d %d", &a, &b);
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    vector<int> order(n);
    for (auto& x: order) scanf("%d", &x);

    DSU dsu(n);
    vector<bool> open(n + 1, false), res(n);
    int components = 0;
    for (int i = n - 1; i >= 0; i--){
        int x = order[i];
        open[x] = true, components++;
        for (int y: adj[x]){
            if (open[y]) components -= dsu.connect(x, y);
        }
        res[i] = components == 1;
    }

    for (int i = 0; i < n; i++) puts(res[i] ? "YES" : "NO");
    return 0;
}
