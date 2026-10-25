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
    int n, m;

    while (scanf("%d %d", &n, &m) == 2){
        DSU dsu(n + m);
        vector<int> node(n + 1), cnt(n + m + 1, 1);
        vector<long long> sum(n + m + 1);
        for (int i = 1; i <= n; i++) node[i] = i, sum[i] = i;

        auto join = [&](int a, int b){
            int c = dsu.find(a), d = dsu.find(b);
            if (!dsu.connect(c, d)) return;
            int r = dsu.find(c);
            cnt[r] = cnt[c] + cnt[d], sum[r] = sum[c] + sum[d];
        };

        for (int fresh = n + 1; m--; ){
            int op, p, q;
            scanf("%d %d", &op, &p);
            if (op == 3){
                int r = dsu.find(node[p]);
                printf("%d %lld\n", cnt[r], sum[r]);
                continue;
            }

            scanf("%d", &q);
            if (op == 1) join(node[p], node[q]);
            else if (!dsu.same(node[p], node[q])){
                int r = dsu.find(node[p]);
                cnt[r]--, sum[r] -= p;
                node[p] = fresh++;
                sum[node[p]] = p;
                join(node[p], node[q]);
            }
        }
    }
    return 0;
}
