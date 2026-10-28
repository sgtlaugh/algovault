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
    freopen("skilevel.in", "r", stdin);
    freopen("skilevel.out", "w", stdout);

    int r, c, t;
    scanf("%d %d %d", &r, &c, &t);

    vector<int> h(r * c), starts(r * c);
    for (auto& x: h) scanf("%d", &x);
    for (auto& x: starts) scanf("%d", &x);

    vector<array<int, 3>> edges;
    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            int id = i * c + j;
            if (j + 1 < c) edges.push_back({abs(h[id] - h[id + 1]), id, id + 1});
            if (i + 1 < r) edges.push_back({abs(h[id] - h[id + c]), id, id + c});
        }
    }
    sort(edges.begin(), edges.end());

    DSU dsu(r * c);
    long long res = 0;
    for (auto& e: edges){
        int a = dsu.find(e[1]), b = dsu.find(e[2]);
        if (a == b) continue;

        int merged = dsu.sz[a] + dsu.sz[b], total = starts[a] + starts[b];
        if (merged >= t){
            if (dsu.sz[a] < t) res += (long long)e[0] * starts[a];
            if (dsu.sz[b] < t) res += (long long)e[0] * starts[b];
        }
        dsu.connect(a, b);
        starts[dsu.find(a)] = total;
    }

    printf("%lld\n", res);
    return 0;
}
