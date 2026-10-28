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
    int n, m, k;

    scanf("%d %d %d", &n, &m, &k);

    vector<long long> keys(m);
    for (auto& key: keys){
        int a, b;
        scanf("%d %d", &a, &b);
        key = (long long)min(a, b) * (n + 1) + max(a, b);
    }
    sort(keys.begin(), keys.end());

    vector<int> order(k);
    vector<bool> broken(m, false);
    for (auto& id: order){
        int a, b;
        scanf("%d %d", &a, &b);
        long long key = (long long)min(a, b) * (n + 1) + max(a, b);
        id = lower_bound(keys.begin(), keys.end(), key) - keys.begin();
        broken[id] = true;
    }

    DSU dsu(n);
    int components = n;
    for (int i = 0; i < m; i++){
        if (!broken[i]) components -= dsu.connect(keys[i] / (n + 1), keys[i] % (n + 1));
    }

    vector<int> res(k);
    for (int i = k - 1; i >= 0; i--){
        res[i] = components;
        components -= dsu.connect(keys[order[i]] / (n + 1), keys[order[i]] % (n + 1));
    }

    for (int i = 0; i < k; i++) printf("%d%c", res[i], i + 1 == k ? '\n' : ' ');
    return 0;
}
