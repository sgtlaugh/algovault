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
    char line[100];
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);
        fgets(line, sizeof(line), stdin);

        DSU dsu(n);
        int successful = 0, unsuccessful = 0;
        /// A blank line or the end of input closes the log of a case
        while (fgets(line, sizeof(line), stdin)){
            char type;
            int a, b;
            if (sscanf(line, " %c %d %d", &type, &a, &b) != 3) break;
            if (type == 'c') dsu.connect(a, b);
            else if (dsu.same(a, b)) successful++;
            else unsuccessful++;
        }
        printf("%d,%d\n", successful, unsuccessful);
        if (tc != t) puts("");
    }
    return 0;
}
