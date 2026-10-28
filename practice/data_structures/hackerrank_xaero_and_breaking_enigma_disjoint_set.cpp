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

const int MAXN = 100010;

int lc[MAXN], rc[MAXN], sub[MAXN];
unsigned pri[MAXN];
bool flip[MAXN];

void push(int t){
    if (t && flip[t]){
        swap(lc[t], rc[t]);
        flip[lc[t]] ^= 1, flip[rc[t]] ^= 1, flip[t] = 0;
    }
}

void pull(int t){
    sub[t] = 1 + sub[lc[t]] + sub[rc[t]];
}

void split(int t, int k, int& a, int& b){
    if (!t){
        a = b = 0;
        return;
    }
    push(t);
    if (sub[lc[t]] < k) split(rc[t], k - sub[lc[t]] - 1, rc[t], b), a = t;
    else split(lc[t], k, a, lc[t]), b = t;
    pull(t);
}

int merge(int a, int b){
    if (!a || !b) return a ? a : b;
    if (pri[a] > pri[b]){
        push(a);
        rc[a] = merge(rc[a], b);
        pull(a);
        return a;
    }
    push(b);
    lc[b] = merge(a, lc[b]);
    pull(b);
    return b;
}

void collect(int t, vector<int>& out){
    if (!t) return;
    push(t);
    collect(lc[t], out);
    out.push_back(t - 1);
    collect(rc[t], out);
}

int main(){
    int n, m;
    long long k;
    static char s[MAXN];

    scanf("%d %d %lld %s", &n, &m, &k, s);

    mt19937 rng(20020523);
    int root = 0;
    for (int i = 1; i <= n; i++){
        pri[i] = rng(), sub[i] = 1;
        root = merge(root, i);
    }

    while (m--){
        int l, r, a, b, c;
        scanf("%d %d", &l, &r);
        split(root, r, b, c);
        split(b, l - 1, a, b);
        flip[b] ^= 1;
        root = merge(merge(a, b), c);
    }

    vector<int> pos;
    collect(root, pos);

    DSU dsu(n);
    for (int i = 0; i < n; i++) dsu.connect(i, pos[i]);

    vector<char> letter(n, '?');
    for (int i = 0; i < n; i++){
        if (s[i] == '?') continue;
        char& c = letter[dsu.find(i)];
        if (c != '?' && c != s[i]){
            puts("Bad Luck Allen");
            return 0;
        }
        c = s[i];
    }

    vector<int> free_roots;
    for (int i = 0; i < n; i++){
        int r = dsu.find(i);
        if (letter[r] == '?') free_roots.push_back(r), letter[r] = 0;
    }

    k--;
    for (int i = (int)free_roots.size() - 1; i >= 0; i--){
        letter[free_roots[i]] = 'a' + k % 26;
        k /= 26;
    }
    if (k){
        puts("Bad Luck Allen");
        return 0;
    }

    for (int i = 0; i < n; i++) s[i] = letter[dsu.find(i)];
    puts(s);
    return 0;
}
