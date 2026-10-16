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
    static char s[200005];
    scanf("%d %d %s", &n, &m, s + 1);
    vector<vector<int>> adj(n + 1);
    for (int i = 0; i < m; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u].push_back(v), adj[v].push_back(u);
    }

    /// removing a node with s = 1 keeps its neighbours connected, so it acts as a permanent relay;
    /// before timestep t the graph is every s = 1 node plus s = 0 nodes >= t, and only nodes >= t count
    DSU dsu(n);
    vector<char> present(n + 1, 0);
    vector<long long> weight(n + 1, 0);
    long long pairs = 0;
    auto join = [&](int a, int b){
        a = dsu.find(a), b = dsu.find(b);
        if (!dsu.connect(a, b)) return;
        int root = dsu.find(a);
        pairs -= weight[a] * (weight[a] - 1) / 2 + weight[b] * (weight[b] - 1) / 2;
        weight[root] = weight[a] + weight[b];
        pairs += weight[root] * (weight[root] - 1) / 2;
    };
    auto add_edges = [&](int v){
        present[v] = 1;
        for (int u : adj[v]){
            if (present[u]) join(u, v);
        }
    };
    for (int v = 1; v <= n; v++){
        if (s[v] == '1') add_edges(v);
    }

    vector<long long> res(n + 1);
    for (int t = n; t >= 1; t--){
        if (s[t] == '0') add_edges(t);
        int root = dsu.find(t);
        pairs += weight[root]++;
        res[t] = pairs;
    }
    for (int t = 1; t <= n; t++) printf("%lld\n", res[t]);
    return 0;
}
