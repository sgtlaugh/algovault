/***
 *
 * Link-Cut Tree
 * Dynamic rooted forest with link, cut, re-rooting, LCA and path aggregates
 *
 * Complexity: amortized O(log n) per operation
 *
 * LinkCutTree<T> lct(n): nodes 0..n - 1 and no edges, values start at T() and paths fold with +
 * LinkCutTree<T, Combine> lct(n, identity, combine): values start at identity and paths fold with combine
 * lct.link(u, v): adds edge u-v with u as a child of v, false and no change if u and v are already connected
 * lct.cut(u, v): removes edge u-v, false if there is no such edge
 * lct.connected(u, v): whether u and v are in the same tree
 * lct.find_root(u): root of u's tree
 * lct.parent(u): parent of u under the current root, -1 for the root
 * lct.make_root(u): re-roots u's tree at u
 * lct.lca(u, v): lowest common ancestor under the current root, -1 if u and v are not connected
 * lct.set(u, value): sets the value of u
 * lct.query(u, v): combine of the values on the path u..v, both ends included, u and v must be connected
 *
 * combine must be associative and commutative with identity as its neutral element, since re-rooting reverses paths
 * link re-roots u's tree at u and v's tree keeps its root, so link(child, parent) builds a rooted forest for lca
 * After cut the piece holding the old root keeps it, the other piece is rooted at its endpoint of the edge
 * query leaves the root where it was
 * Everything is iterative, deep paths are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Combine = plus<T>>
struct LinkCutTree{
    struct Node{
        /// par of an auxiliary splay tree's root is the path-parent pointer into the splay tree above it
        int ch[2], par;
        bool flip;
        T value, agg;
    };

    T identity;
    Combine combine;
    vector<Node> nodes;

    LinkCutTree(int n, T identity = T(), Combine combine = Combine())
        : identity(identity), combine(combine), nodes(n, Node{{-1, -1}, -1, false, identity, identity}) {}

    int parent(int u){
        access(u);
        int x = nodes[u].ch[0];
        if (x == -1) return -1;

        for (push(x); nodes[x].ch[1] != -1; push(x)) x = nodes[x].ch[1];
        splay(x);
        return x;
    }

    int find_root(int u){
        access(u);
        int x = u;
        for (push(x); nodes[x].ch[0] != -1; push(x)) x = nodes[x].ch[0];

        splay(x);
        return x;
    }

    bool connected(int u, int v){
        return find_root(u) == find_root(v);
    }

    void make_root(int u){
        access(u);
        nodes[u].flip ^= 1;
    }

    bool link(int u, int v){
        if (connected(u, v)) return false;
        make_root(u);
        nodes[u].par = v;
        return true;
    }

    bool cut(int u, int v){
        if (parent(u) != v) swap(u, v);
        if (parent(u) != v) return false;

        access(u);
        nodes[nodes[u].ch[0]].par = -1;
        nodes[u].ch[0] = -1;
        pull(u);
        return true;
    }

    int lca(int u, int v){
        if (!connected(u, v)) return -1;
        access(u);
        return access(v);
    }

    void set(int u, T value){
        splay(u);
        nodes[u].value = value;
        pull(u);
    }

    T query(int u, int v){
        assert(connected(u, v));
        int root = find_root(u);

        make_root(u);
        access(v);
        T res = nodes[v].agg;

        make_root(root);
        return res;
    }

    bool is_root(int x) const{
        int p = nodes[x].par;
        return p == -1 || (nodes[p].ch[0] != x && nodes[p].ch[1] != x);
    }

    T agg_of(int x) const{
        return x == -1 ? identity : nodes[x].agg;
    }

    void push(int x){
        if (!nodes[x].flip) return;
        swap(nodes[x].ch[0], nodes[x].ch[1]);
        for (int c : nodes[x].ch){
            if (c != -1) nodes[c].flip ^= 1;
        }
        nodes[x].flip = false;
    }

    void pull(int x){
        nodes[x].agg = combine(combine(agg_of(nodes[x].ch[0]), nodes[x].value), agg_of(nodes[x].ch[1]));
    }

    void rotate(int x){
        int p = nodes[x].par, g = nodes[p].par, d = nodes[p].ch[1] == x, child = nodes[x].ch[!d];

        if (!is_root(p)) nodes[g].ch[nodes[g].ch[1] == p] = x;
        nodes[x].par = g;
        nodes[p].ch[d] = child;
        if (child != -1) nodes[child].par = p;
        nodes[x].ch[!d] = p;
        nodes[p].par = x;

        pull(p);
        pull(x);
    }

    /// Pushing only g, p, x is enough: a flip pending higher up reverses the whole subtree and commutes with rotations below
    void splay(int x){
        push(x);
        while (!is_root(x)){
            int p = nodes[x].par, g = nodes[p].par;
            if (!is_root(p)) push(g);
            push(p);
            push(x);
            if (!is_root(p)) rotate((nodes[g].ch[1] == p) == (nodes[p].ch[1] == x) ? p : x);
            rotate(x);
        }
    }

    /// Makes root..x the preferred path, x ends as the root of its splay tree; returns the last node joined
    int access(int x){
        int last = -1;
        for (int y = x; y != -1; y = nodes[y].par){
            splay(y);
            nodes[y].ch[1] = last;
            pull(y);
            last = y;
        }

        splay(x);
        return last;
    }
};

int main(){
    int n, m;
    scanf("%d %d", &n, &m);

    /// Edges become nodes n..n + m - 1 holding weight * 2^20 + index, distinct weights make the path max name one edge
    auto heaviest = [](long long a, long long b){ return max(a, b); };
    LinkCutTree<long long, decltype(heaviest)> lct(n + m, -1, heaviest);
    vector<int> eu(m), ev(m);
    for (int i = 0; i < m; i++){
        long long w;
        scanf("%d %d %lld", &eu[i], &ev[i], &w);
        int u = eu[i], v = ev[i], e = n + i;

        int res = -1;
        if (lct.connected(u, v)){
            long long top = lct.query(u, v);
            if (top < (w << 20)){
                res = i;
            }
            else{
                res = top & ((1 << 20) - 1);
                lct.cut(n + res, eu[res]);
                lct.cut(n + res, ev[res]);
            }
        }

        if (res != i){
            lct.set(e, w << 20 | i);
            lct.link(e, u);
            lct.link(v, e);
        }
        printf("%d%c", res, i + 1 == m ? '\n' : ' ');
    }
    return 0;
}
