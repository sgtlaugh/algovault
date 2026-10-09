/***
 *
 * Link-Cut Tree
 * Dynamic rooted forest on nodes 0..n-1 with link, cut, re-rooting, LCA and path aggregates
 *
 * Complexity: amortized O(log n) per operation, O(n) memory
 *
 * LinkCutTree<T, Combine> lct(n, identity, combine): n isolated nodes, each starting with value identity
 *   link(u, v): adds edge u-v with u as a child of v (u's tree is re-rooted at u, v's tree keeps its root)
 *               false and no change if u and v are already connected
 *   cut(u, v): removes edge u-v, false if there is no such edge
 *              the piece holding the old root keeps it, the other piece is rooted at its endpoint of the edge
 *   connected(u, v): whether u and v are in the same tree
 *   find_root(u): root of u's tree
 *   parent(u): parent of u under the current root, -1 for the root
 *   make_root(u): re-roots u's tree at u
 *   lca(u, v): lowest common ancestor under the current root, -1 if not connected
 *   set(u, value): point update
 *   query(u, v): combine of the values on the path u..v (both ends included), u and v must be connected
 *                the tree's root is restored afterwards
 *
 * combine must be associative and commutative with identity as its neutral element, since re-rooting reverses paths
 *   sum: LinkCutTree<long long> lct(n);
 *   max: auto mx = [](long long a, long long b){ return max(a, b); };
 *        LinkCutTree<long long, decltype(mx)> lct(n, LLONG_MIN, mx);
 *
 * link(child, parent) builds a rooted forest, so lca works without any make_root call
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
};

int main(){
    LinkCutTree<long long> sum(7);
    auto mx = [](long long a, long long b){ return max(a, b); };
    LinkCutTree<long long, decltype(mx)> best(7, LLONG_MIN, mx);

    vector<long long> values = {5, 3, 8, 1, 7, 2, 4};
    for (int i = 0; i < 7; i++) sum.set(i, values[i]), best.set(i, values[i]);
    for (auto [u, v] : vector<pair<int, int>>{{1, 0}, {2, 0}, {3, 1}, {4, 1}, {5, 2}}){
        assert(sum.link(u, v) && best.link(u, v));
    }

    assert(sum.connected(3, 5) && !sum.connected(3, 6) && sum.connected(6, 6));
    assert(sum.find_root(4) == 0 && sum.find_root(6) == 6);
    assert(sum.parent(3) == 1 && sum.parent(5) == 2 && sum.parent(0) == -1 && sum.parent(6) == -1);
    assert(sum.lca(3, 4) == 1 && sum.lca(3, 5) == 0 && sum.lca(4, 1) == 1 && sum.lca(2, 2) == 2 && sum.lca(3, 6) == -1);
    assert(sum.query(3, 5) == 19 && best.query(3, 5) == 8);
    assert(sum.query(4, 4) == 7 && best.query(4, 3) == 7);
    assert(sum.find_root(3) == 0);
    assert(!sum.link(3, 5) && !sum.link(0, 0));

    sum.make_root(2);
    assert(sum.find_root(3) == 2 && sum.lca(3, 4) == 1 && sum.lca(0, 5) == 2 && sum.lca(3, 5) == 2);
    assert(sum.parent(0) == 2 && sum.parent(2) == -1 && sum.parent(5) == 2 && sum.parent(1) == 0);

    assert(sum.cut(0, 1) && !sum.cut(0, 1) && !sum.cut(3, 4) && !sum.cut(3, 6));
    assert(!sum.connected(3, 5) && sum.find_root(4) == 1 && sum.find_root(0) == 2);
    sum.set(1, 100);
    assert(sum.query(3, 4) == 108);

    assert(sum.link(1, 5));
    assert(sum.find_root(3) == 2 && sum.lca(3, 0) == 2 && sum.lca(4, 5) == 5);
    assert(sum.query(4, 0) == 122);

    LinkCutTree<long long, decltype(mx)> negative(3, LLONG_MIN, mx);
    negative.set(0, -5), negative.set(1, -2), negative.set(2, -9);
    assert(negative.link(0, 1) && negative.link(2, 1));
    assert(negative.query(0, 2) == -2 && negative.query(0, 0) == -5 && negative.query(2, 2) == -9);
    assert(negative.cut(1, 2) && negative.query(0, 1) == -2 && negative.find_root(2) == 2);

    return 0;
}
