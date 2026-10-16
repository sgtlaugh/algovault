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

const int MAXV = 100000;

int main(){
    int t;
    scanf("%d", &t);

    for (int tc = 1; tc <= t; tc++){
        int n, q;
        scanf("%d %d", &n, &q);

        /// Positions sharing a value form one set; holder[x] is a node of x's set, value[root] its value
        DSU dsu(n);
        vector<int> holder(MAXV + 1, 0), value(n + 1);
        for (int i = 1; i <= n; i++){
            scanf("%d", &value[i]);
            if (holder[value[i]]) dsu.connect(holder[value[i]], i);
            holder[value[i]] = i;
        }
        for (int x = 1; x <= MAXV; x++){
            if (holder[x]) value[dsu.find(holder[x])] = x;
        }

        printf("Case %d:\n", tc);
        while (q--){
            int type, x, y;
            scanf("%d %d", &type, &x);
            if (type == 2){
                printf("%d\n", value[dsu.find(x)]);
                continue;
            }

            scanf("%d", &y);
            if (x == y || !holder[x]) continue;
            if (holder[y]) dsu.connect(holder[x], holder[y]);
            else holder[y] = holder[x];
            holder[x] = 0;
            value[dsu.find(holder[y])] = y;
        }
    }
    return 0;
}
