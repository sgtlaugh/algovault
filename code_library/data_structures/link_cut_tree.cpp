/***
 *
 * Link-Cut Tree
 * Dynamic rooted forest with link, cut, re-rooting, LCA and path aggregates
 *
 * Complexity: amortized O(log n) per operation
 *
 * LinkCutTree<T = long long, Ops = PathSumAdd<T>, L = T> lct(n): nodes 0..n - 1 and no edges, values start at T()
 * lct.cut(u, v): removes edge u-v, false if there is no such edge
 * lct.link(u, v): adds edge u-v with u as a child of v, false and no change if u and v are already connected
 * lct.make_root(u): re-roots u's tree at u
 * lct.set(u, value): sets the value of u
 * lct.update(u, v, lz): applies the lazy update lz to every value on the path u..v, adds lz by default
 * lct.connected(u, v): whether u and v are in the same tree
 * lct.find_root(u): root of u's tree
 * lct.lca(u, v): lowest common ancestor under the current root, -1 if u and v are not connected
 * lct.parent(u): parent of u under the current root, -1 for the root
 * lct.query(u, v): combine of the values on the path from u to v in that order, ends included, their sum by default
 *
 * Notes:
 *   - update and query need u and v connected, and both leave the root where it was
 *   - Ops holds four static functions, like SegmentTree's, L is the lazy type when it differs from T:
 *       - combine(a, b): fold of two adjacent path pieces, a nearer u
 *       - apply(val, lz, len): value of a piece of len nodes after lz, len = 1 for a single value
 *       - compose(old, new): one lazy meaning old then new, also called on a stale old whose result is dropped
 *       - reversed(a): the aggregate of the piece read backwards, a itself when combine is commutative
 *   - link re-roots u's tree at u and v's tree keeps its root, so link(child, parent) builds a rooted forest for lca
 *   - After cut the piece holding the old root keeps it, the other piece is rooted at its endpoint of the edge
 *   - n < 2^30, and everything is iterative, so deep paths are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Path add and path sum, the default Ops of LinkCutTree
template <typename T>
struct PathSumAdd{
    static T combine(T a, T b){ return a + b; }
    static T apply(T val, T lz, int len){ return val + lz * len; }
    static T compose(T old_lz, T new_lz){ return old_lz + new_lz; }
    static T reversed(T a){ return a; }
};

template <typename T = long long, typename Ops = PathSumAdd<T>, typename L = T>
struct LinkCutTree{
    struct Node{
        /// par of an auxiliary splay tree's root is the path-parent pointer into the splay tree above it
        int ch[2], par;
        unsigned meta;
        T value, agg;
    };

    /// meta packs the splay subtree size with the two pending flags, which keeps a node at 32 bytes for T = long long
    static const unsigned FLIP = 1u << 30, LAZY = 1u << 31, SIZE = FLIP - 1;
    vector<Node> nodes;
    vector<L> lazy;

    LinkCutTree(int n) : nodes(n, Node{{-1, -1}, -1, 1, T(), T()}), lazy(n){
        assert(0 <= n && n < (1 << 30));
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

    bool link(int u, int v){
        if (connected(u, v)) return false;
        make_root(u);
        nodes[u].par = v;
        return true;
    }

    void make_root(int u){
        access(u);
        apply_flip(u);
    }

    void set(int u, const T& value){
        splay(u);
        nodes[u].value = value;
        pull(u);
    }

    void update(int u, int v, const L& lz){
        int w = meet(u, v);
        assert(w != -1);

        splay(w);
        apply_lazy(nodes[w].ch[1], lz);
        nodes[w].value = Ops::apply(nodes[w].value, lz, 1);
        pull(w);
        if (w != u) apply_lazy(u, lz);
    }

    bool connected(int u, int v){
        return meet(u, v) != -1;
    }

    int find_root(int u){
        access(u);
        int x = u;
        for (push(x); nodes[x].ch[0] != -1; push(x)) x = nodes[x].ch[0];

        splay(x);
        return x;
    }

    int lca(int u, int v){
        return meet(u, v);
    }

    int parent(int u){
        access(u);
        int x = nodes[u].ch[0];
        if (x == -1) return -1;

        for (push(x); nodes[x].ch[1] != -1; push(x)) x = nodes[x].ch[1];
        splay(x);
        return x;
    }

    T query(int u, int v){
        int w = meet(u, v);
        assert(w != -1);

        splay(w);
        T res = nodes[w].value;
        if (nodes[w].ch[1] != -1) res = Ops::combine(res, nodes[nodes[w].ch[1]].agg);
        if (w != u) res = Ops::combine(Ops::reversed(nodes[u].agg), res);
        return res;
    }

    /// LCA of u and v, -1 if not connected; when it is not u, u ends as the root of the splay tree of the path below it
    int meet(int u, int v){
        access(u);
        int w = access(v);
        if (w == u) return w;

        splay(u);
        return nodes[u].par == -1 ? -1 : w;
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

    void pull(int x){
        Node& node = nodes[x];
        T agg = node.value;
        unsigned size = 1;
        if (node.ch[0] != -1) agg = Ops::combine(nodes[node.ch[0]].agg, agg), size += nodes[node.ch[0]].meta & SIZE;
        if (node.ch[1] != -1) agg = Ops::combine(agg, nodes[node.ch[1]].agg), size += nodes[node.ch[1]].meta & SIZE;
        node.agg = agg, node.meta = (node.meta & ~SIZE) | size;
    }

    /// Pushing only g, p, x is enough: an update pending higher up covers the whole subtree, which rotations keep
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

    void apply_lazy(int x, const L& lz){
        if (x == -1) return;
        Node& node = nodes[x];
        node.value = Ops::apply(node.value, lz, 1);
        node.agg = Ops::apply(node.agg, lz, node.meta & SIZE);
        /// compose also runs on a stale lazy and the result is dropped: a select, not a branch, like ImplicitTreap
        L composed = Ops::compose(lazy[x], lz);
        lazy[x] = (node.meta & LAZY) ? composed : lz;
        node.meta |= LAZY;
    }

    /// The children swap at once and agg is read backwards, so a flipped node's own fields are already correct
    void apply_flip(int x){
        if (x == -1) return;
        Node& node = nodes[x];
        swap(node.ch[0], node.ch[1]);
        node.agg = Ops::reversed(node.agg);
        node.meta ^= FLIP;
    }

    void push(int x){
        Node& node = nodes[x];
        if (!(node.meta & (FLIP | LAZY))) return;

        int l = node.ch[0], r = node.ch[1];
        if (node.meta & LAZY){
            apply_lazy(l, lazy[x]);
            apply_lazy(r, lazy[x]);
            node.meta &= ~LAZY;
        }
        if (node.meta & FLIP){
            apply_flip(l);
            apply_flip(r);
            node.meta &= ~FLIP;
        }
    }

    bool is_root(int x) const{
        int p = nodes[x].par;
        return p == -1 || (nodes[p].ch[0] != x && nodes[p].ch[1] != x);
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
};

int main(){
    /// Path sums on a forest rooted at 0 and 6
    /***
     *        0           6
     *      /   \
     *     1     2
     *    / \     \
     *   3   4     5
    ***/
    LinkCutTree<long long> sum(7);
    vector<long long> values = {5, 3, 8, 1, 7, 2, 4};
    for (int i = 0; i < 7; i++) sum.set(i, values[i]);
    assert(sum.link(1, 0));
    assert(sum.link(2, 0));
    assert(sum.link(3, 1));
    assert(sum.link(4, 1));
    assert(sum.link(5, 2));
    assert(!sum.link(3, 5));

    assert(sum.find_root(4) == 0);
    assert(!sum.connected(3, 6));
    assert(sum.parent(3) == 1);
    assert(sum.lca(3, 5) == 0);
    assert(sum.query(3, 5) == 19);              /// path 3 1 0 2 5: 1 + 3 + 5 + 8 + 2

    sum.make_root(2);
    assert(sum.lca(3, 4) == 1);
    assert(sum.parent(0) == 2);

    assert(sum.cut(0, 1));
    assert(sum.find_root(4) == 1);              /// the piece without the root is rooted at the cut endpoint

    /// Path maximum with path add through an Ops struct, max reads the same both ways
    struct MaxAdd{
        static long long combine(long long a, long long b){ return max(a, b); }
        static long long apply(long long val, long long lz, int){ return val + lz; }
        static long long compose(long long old_lz, long long new_lz){ return old_lz + new_lz; }
        static long long reversed(long long a){ return a; }
    };

    LinkCutTree<long long, MaxAdd> best(3);
    best.set(0, -5);
    best.set(1, -2);
    best.set(2, -9);
    assert(best.link(0, 1));
    assert(best.link(2, 1));
    assert(best.query(0, 2) == -2);
    best.update(0, 1, 10);                      /// -5 and -2 become 5 and 8
    assert(best.query(0, 2) == 8);
    assert(best.query(2, 2) == -9);

    return 0;
}
