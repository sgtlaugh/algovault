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
    freopen("fcolor.in", "r", stdin);
    freopen("fcolor.out", "w", stdout);

    int n, m;
    scanf("%d %d", &n, &m);
    vector<vector<int>> admirers(n + 1);
    for (int i = 0; i < m; i++){
        int a, b;
        scanf("%d %d", &a, &b);
        admirers[a].push_back(b);
    }

    /// all admirers of one color class share a color, merging them merges their admirer lists, repeat until stable
    DSU dsu(n);
    queue<int> pending;
    for (int i = 1; i <= n; i++){
        if (admirers[i].size() > 1) pending.push(i);
    }
    while (!pending.empty()){
        int r = pending.front();
        pending.pop();
        while (admirers[r].size() > 1){
            int a = dsu.find(admirers[r].back());
            admirers[r].pop_back();
            int b = dsu.find(admirers[r].back());
            if (!dsu.connect(a, b)) continue;

            int root = dsu.find(a), other = root == a ? b : a;
            if (admirers[root].size() < admirers[other].size()) swap(admirers[root], admirers[other]);
            admirers[root].insert(admirers[root].end(), admirers[other].begin(), admirers[other].end());
            vector<int>().swap(admirers[other]);
            pending.push(root);
        }
    }

    vector<int> color(n + 1, 0);
    int colors = 0;
    for (int i = 1; i <= n; i++){
        int root = dsu.find(i);
        if (!color[root]) color[root] = ++colors;
        printf("%d\n", color[root]);
    }
    return 0;
}
