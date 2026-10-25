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
    freopen("wormsort.in", "r", stdin);
    freopen("wormsort.out", "w", stdout);

    int n, m;
    scanf("%d %d", &n, &m);
    vector<int> p(n + 1);
    for (int i = 1; i <= n; i++) scanf("%d", &p[i]);
    vector<array<int, 3>> edges(m);
    for (auto& e : edges) scanf("%d %d %d", &e[1], &e[2], &e[0]);
    sort(edges.rbegin(), edges.rend());

    DSU dsu(n);
    int ptr = 1;
    auto advance = [&](){ while (ptr <= n && dsu.same(ptr, p[ptr])) ptr++; };
    advance();
    if (ptr > n){
        puts("-1");
        return 0;
    }
    for (auto& e : edges){
        dsu.connect(e[1], e[2]);
        advance();
        if (ptr > n){
            printf("%d\n", e[0]);
            return 0;
        }
    }
    return 0;
}
