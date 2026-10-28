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
    freopen("valleys.in", "r", stdin);
    freopen("valleys.out", "w", stdout);

    int n;
    scanf("%d", &n);
    vector<int> h(n * n), order(n * n);
    for (int& x : h) scanf("%d", &x);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){ return h[a] < h[b]; });

    DSU dsu(n * n);
    vector<char> on(n * n, 0);
    vector<int> chi(n * n, 0);
    const int dr[] = {-1, 1, 0, 0}, dc[] = {0, 0, -1, 1};
    auto added = [&](int r, int c){ return r >= 0 && r < n && c >= 0 && c < n && on[r * n + c]; };

    long long res = 0;
    for (int cell : order){
        int r = cell / n, c = cell % n, edges = 0, faces = 0;
        on[cell] = 1, chi[cell] = 1;
        for (int d = 0; d < 4; d++){
            if (!added(r + dr[d], c + dc[d])) continue;
            int nb = (r + dr[d]) * n + c + dc[d], ra = dsu.find(cell), rb = dsu.find(nb);
            edges++;
            if (dsu.connect(ra, rb)) chi[dsu.find(cell)] = chi[ra] + chi[rb];
        }
        for (int a = -1; a <= 1; a += 2){
            for (int b = -1; b <= 1; b += 2){
                if (added(r + a, c) && added(r, c + b) && added(r + a, c + b)) faces++;
            }
        }
        int root = dsu.find(cell);
        chi[root] += faces - edges;
        if (chi[root] == 1) res += dsu.size(root);
    }

    printf("%lld\n", res);
    return 0;
}
