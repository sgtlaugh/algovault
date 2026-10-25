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

const int MAXH = 7000;

int main(){
    int n, m;

    scanf("%d %d", &n, &m);

    vector<array<int, 3>> edges(m);
    for (auto& e: edges) scanf("%d %d %d", &e[1], &e[2], &e[0]);
    sort(edges.rbegin(), edges.rend());

    vector<int> cnt(n + 1, 0);
    set<int> sizes = {1};
    cnt[1] = n;

    auto remove_size = [&](int s){
        if (--cnt[s] == 0) sizes.erase(s);
    };

    DSU dsu(n - 1);
    double res = 1e18;
    for (auto& e: edges){
        int a = dsu.find(e[1]), b = dsu.find(e[2]);
        if (a == b) continue;

        bitset<MAXH + 1> reach;
        reach[0] = 1;
        for (int s: sizes){
            int c = cnt[s];
            for (int take = 1; c > 0; take <<= 1){
                int piece = min(take, c);
                c -= piece;
                if ((long long)piece * s <= MAXH) reach |= reach << (piece * s);
            }
        }
        int half = n / 2;
        while (!reach[half]) half--;
        res = min(res, (double)e[0] / half);

        int sa = dsu.sz[a], sb = dsu.sz[b];
        remove_size(sa), remove_size(sb);
        dsu.connect(a, b);
        if (cnt[sa + sb]++ == 0) sizes.insert(sa + sb);
    }

    printf("%.12f\n", res);
    return 0;
}
