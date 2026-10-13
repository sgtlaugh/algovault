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
    int n;
    scanf("%d", &n);

    DSU dsu(n);
    vector<pair<int, int>> extra;
    for (int i = 1; i < n; i++){
        int a, b;
        scanf("%d %d", &a, &b);
        /// A road inside an already connected component closes a cycle, so dropping it disconnects nothing
        if (!dsu.connect(a, b)) extra.push_back({a, b});
    }

    vector<int> roots;
    for (int i = 1; i <= n; i++){
        if (dsu.find(i) == i) roots.push_back(i);
    }

    /// n - 1 roads over c components leave exactly c - 1 cycle roads, one per missing link between components
    printf("%d\n", (int)extra.size());
    for (int i = 0; i < (int)extra.size(); i++){
        printf("%d %d %d %d\n", extra[i].first, extra[i].second, roots[0], roots[i + 1]);
    }
    return 0;
}
