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

    vector<int> u(m + 1), v(m + 1);
    for (int i = 1; i <= m; i++) scanf("%d %d", &u[i], &v[i]);

    int q;
    scanf("%d", &q);
    vector<int> order(q);
    vector<bool> torn(m + 1, false);
    for (auto& id: order){
        scanf("%d", &id);
        torn[id] = true;
    }

    /// Tearing only splits pieces, so replay it backwards as unions
    DSU dsu(n);
    int pieces = n;
    for (int i = 1; i <= m; i++){
        if (!torn[i]) pieces -= dsu.connect(u[i], v[i]);
    }

    vector<int> res(q);
    for (int i = q - 1; i >= 0; i--){
        res[i] = pieces;
        pieces -= dsu.connect(u[order[i]], v[order[i]]);
    }

    for (int i = 0; i < q; i++) printf("%d%c", res[i], i + 1 == q ? '\n' : ' ');
    return 0;
}
