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
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        double d;
        scanf("%d %lf", &n, &d);
        /// Two decimals at most, so hundredths make every distance comparison exact
        long long r = llround(d * 100);
        vector<long long> x(n), y(n);
        for (int i = 0; i < n; i++){
            double a, b;
            scanf("%lf %lf", &a, &b);
            x[i] = llround(a * 100), y[i] = llround(b * 100);
        }

        DSU dsu(n);
        int res = n;
        for (int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                long long dx = x[i] - x[j], dy = y[i] - y[j];
                if (dx * dx + dy * dy <= r * r && dsu.connect(i, j)) res--;
            }
        }
        printf("Case %d: %d\n", tc, res);
    }
    return 0;
}
