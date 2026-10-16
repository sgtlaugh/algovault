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
    freopen("moocast.in", "r", stdin);
    freopen("moocast.out", "w", stdout);

    int n;
    scanf("%d", &n);
    vector<int> x(n), y(n);
    for (int i = 0; i < n; i++) scanf("%d %d", &x[i], &y[i]);

    vector<tuple<int, int, int>> edges;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            int dx = x[i] - x[j], dy = y[i] - y[j];
            edges.emplace_back(dx * dx + dy * dy, i, j);
        }
    }
    sort(edges.begin(), edges.end());

    DSU dsu(n);
    int components = n, res = 0;
    for (auto& e : edges){
        if (components == 1) break;
        if (dsu.connect(get<1>(e), get<2>(e))) components--, res = get<0>(e);
    }
    printf("%d\n", res);
    return 0;
}
