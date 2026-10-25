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
    int n, m, k;

    scanf("%d %d %d", &n, &m, &k);

    vector<int> gov(k);
    for (auto& c : gov) scanf("%d", &c);

    DSU dsu(n);
    for (int i = 0; i < m; i++){
        int u, v;
        scanf("%d %d", &u, &v);
        dsu.connect(u, v);
    }

    long long free_nodes = n, largest = 0, res = 0;
    for (int c : gov){
        long long s = dsu.size(c);
        free_nodes -= s;
        largest = max(largest, s);
        res += s * (s - 1) / 2;
    }

    long long merged = largest + free_nodes;
    res += merged * (merged - 1) / 2 - largest * (largest - 1) / 2;
    printf("%lld\n", res - m);
    return 0;
}
