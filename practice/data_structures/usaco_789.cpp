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
    freopen("mootube.in", "r", stdin);
    freopen("mootube.out", "w", stdout);

    int n, q;
    scanf("%d %d", &n, &q);
    vector<array<int, 3>> edges(n - 1), queries(q);
    for (auto& e : edges) scanf("%d %d %d", &e[1], &e[2], &e[0]);
    for (int i = 0; i < q; i++) scanf("%d %d", &queries[i][0], &queries[i][1]), queries[i][2] = i;
    sort(edges.rbegin(), edges.rend());
    sort(queries.rbegin(), queries.rend());

    DSU dsu(n);
    vector<int> res(q);
    size_t j = 0;
    for (auto& qr : queries){
        while (j < edges.size() && edges[j][0] >= qr[0]) dsu.connect(edges[j][1], edges[j][2]), j++;
        res[qr[2]] = dsu.size(qr[1]) - 1;
    }

    for (int x : res) printf("%d\n", x);
    return 0;
}
