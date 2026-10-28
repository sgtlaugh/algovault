/***
 *
 * Euler Tour Tree
 * Dynamic forest with link, cut, subtree sums and subtree range add
 *
 * Complexity: amortized O(log n) per operation
 *
 * EulerTourTree<T> ett(n): nodes 0..n - 1 with value T() and no edges
 * EulerTourTree<T> ett(values): nodes 0..n - 1 with the given values and no edges
 * ett.add(u, x): adds x to the value of u
 * ett.cut(u, v): removes edge u-v, false if there is no such edge
 * ett.link(u, v): adds edge u-v, false and no change if u and v are already connected
 * ett.set(u, x): sets the value of u to x
 * ett.subtree_add(v, p, x): adds x to every value in v's subtree when its tree is rooted at p, edge v-p must exist
 * ett.connected(u, v): whether u and v are in the same tree
 * ett.get(u): the value of u
 * ett.subtree_query(v, p): sum of the values in v's subtree when its tree is rooted at p, edge v-p must exist
 *
 * Notes:
 *   - T needs + and * by an int count, so long long or a modint; sums must fit in T
 *   - There is no fixed root: every subtree call names the parent p, so any neighbor of v can play the parent
 *   - For path queries use link_cut_tree.cpp, which has no subtree aggregates
 *   - Cut edges leave two dead nodes in nodes, so memory grows by two nodes per link
 *   - Everything is iterative, deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct EulerTourTree{
    /// Nodes 0..n - 1 are the vertices, the rest are arcs u->v with value T(); a tree is its tour in one splay tree
    struct Node{
        int ch[2], par, vertices;
        T value, sum, lazy;
    };

    int n;
    vector<Node> nodes;
    /// Edge min(u, v)-max(u, v) maps to arc min->max, and arc max->min is the next node
    unordered_map<long long, int> edges;

    EulerTourTree(int n) : EulerTourTree(vector<T>(n)) {}

    EulerTourTree(const vector<T>& values) : n(values.size()){
        for (const T& x : values) nodes.push_back({{-1, -1}, -1, 1, x, x, T()});
    }

    void add(int u, T x){
        splay(u);
        nodes[u].value += x;
        pull(u);
    }

    bool cut(int u, int v){
        int uv = arc(u, v), vu = arc(v, u);
        if (uv == -1) return false;
        edges.erase(key(u, v));

        /// The two arcs split the cyclic tour into the part between them and two outer parts that close into one tour
        int d = splay_pair(uv, vu);
        int between = nodes[vu].ch[!d], outer_a = nodes[uv].ch[!d], outer_b = nodes[vu].ch[d];
        for (int c : {between, outer_a, outer_b}){
            if (c != -1) nodes[c].par = -1;
        }
        join(outer_a, outer_b);
        return true;
    }

    bool link(int u, int v){
        if (connected(u, v)) return false;
        int a = nodes.size();
        edges[key(u, v)] = a;
        nodes.resize(a + 2, Node{{-1, -1}, -1, 0, T(), T(), T()});

        /// connected left v at the root: its tour A v B becomes A v (v->u) (u's tour from u) (u->v) B
        int after_v = nodes[v].ch[1];
        if (after_v != -1) nodes[after_v].par = -1;
        nodes[v].ch[1] = -1;

        int tail = join(make_first(u), join(arc(u, v), after_v));
        join(v, join(arc(v, u), tail));
        return true;
    }

    void set(int u, T x){
        splay(u);
        nodes[u].value = x;
        pull(u);
    }

    void subtree_add(int v, int p, T x){
        auto [pv, vp] = split_subtree(v, p);
        apply(nodes[vp].ch[0], x);
        if (nodes[pv].ch[0] == vp) apply(nodes[pv].ch[1], x);
        pull(vp);
        pull(pv);
    }

    bool connected(int u, int v){
        splay(u);
        splay(v);
        return u == v || nodes[u].par != -1;
    }

    T get(int u){
        splay(u);
        return nodes[u].value;
    }

    T subtree_query(int v, int p){
        auto [pv, vp] = split_subtree(v, p);
        T res = sum_of(nodes[vp].ch[0]);
        if (nodes[pv].ch[0] == vp) res += sum_of(nodes[pv].ch[1]);
        return res;
    }

    long long key(int u, int v) const{
        return (long long)min(u, v) * n + max(u, v);
    }

    int arc(int u, int v) const{
        auto it = edges.find(key(u, v));
        if (it == edges.end()) return -1;
        return it->second + (u > v);
    }

    /// Splays a to the root and b to a child of a, both in one tour; returns 1 if b comes after a
    int splay_pair(int a, int b){
        splay(a);
        splay(b, a);
        return nodes[a].ch[1] == b;
    }

    /// The tour is cyclic, so rotating it to start at any node keeps it a valid tour; returns the new root
    int make_first(int x){
        splay(x);
        int left = nodes[x].ch[0];
        if (left == -1) return x;

        nodes[left].par = -1;
        nodes[x].ch[0] = -1;
        pull(x);
        return join(x, left);
    }

    /// Joins the tours rooted at a and b, a first; returns the new root
    int join(int a, int b){
        if (a == -1 || b == -1) return a == -1 ? b : a;
        int x = a;
        for (push(x); nodes[x].ch[1] != -1; push(x)) x = nodes[x].ch[1];

        splay(x);
        nodes[x].ch[1] = b;
        nodes[b].par = x;
        pull(x);
        return x;
    }

    /// v's side is the cyclic tour from arc p->v to arc v->p: vp's left child, plus pv's right child if vp is before pv
    pair<int, int> split_subtree(int v, int p){
        int pv = arc(p, v), vp = arc(v, p);
        assert(pv != -1);

        splay_pair(pv, vp);
        return {pv, vp};
    }

    T sum_of(int x) const{
        return x == -1 ? T() : nodes[x].sum;
    }

    /// Pending lazies above g are uniform over the subtree being rotated, so pushing only g, p, x is enough
    void splay(int x, int goal = -1){
        push(x);
        while (nodes[x].par != goal){
            int p = nodes[x].par, g = nodes[p].par;
            if (g != goal) push(g);
            push(p);
            push(x);
            if (g != goal) rotate((nodes[g].ch[1] == p) == (nodes[p].ch[1] == x) ? p : x);
            rotate(x);
        }
        pull(x);
    }

    void push(int x){
        if (nodes[x].lazy == T()) return;
        for (int c : nodes[x].ch) apply(c, nodes[x].lazy);
        nodes[x].lazy = T();
    }

    void apply(int x, T d){
        if (x == -1) return;
        if (x < n) nodes[x].value += d;
        nodes[x].sum += d * nodes[x].vertices;
        nodes[x].lazy += d;
    }

    /// Pulls only p, x stays stale until splay pulls it once at the end
    void rotate(int x){
        int p = nodes[x].par, g = nodes[p].par, d = nodes[p].ch[1] == x, child = nodes[x].ch[!d];

        if (g != -1) nodes[g].ch[nodes[g].ch[1] == p] = x;
        nodes[x].par = g;
        nodes[p].ch[d] = child;
        if (child != -1) nodes[child].par = p;
        nodes[x].ch[!d] = p;
        nodes[p].par = x;

        pull(p);
    }

    void pull(int x){
        nodes[x].vertices = x < n;
        nodes[x].sum = nodes[x].value;
        for (int c : nodes[x].ch){
            if (c == -1) continue;
            nodes[x].vertices += nodes[c].vertices;
            nodes[x].sum += nodes[c].sum;
        }
    }
};

int main(){
    /// Subtree sums on the tree below, asked with respect to different parents
    /***
     *        0
     *      /   \
     *     1     2
     *    / \
     *   3   4
    ***/
    EulerTourTree<long long> ett(vector<long long>{5, 3, 8, 1, 7});
    assert(ett.link(1, 0));
    assert(ett.link(2, 0));
    assert(ett.link(3, 1));
    assert(ett.link(4, 1));
    assert(!ett.link(3, 2));

    assert(ett.subtree_query(1, 0) == 11);          /// 3 + 1 + 7
    assert(ett.subtree_query(0, 1) == 13);          /// rooted at 1, 0's side is 0 and 2: 5 + 8
    assert(ett.subtree_query(3, 1) == 1);

    ett.subtree_add(1, 0, 10);                      /// 1, 3 and 4 gain 10
    assert(ett.get(4) == 17);
    assert(ett.subtree_query(0, 1) == 13);
    ett.add(2, 2);
    assert(ett.subtree_query(0, 1) == 15);

    /// Moving a subtree: 4 hangs from 2 instead of 1
    assert(ett.cut(4, 1));
    assert(!ett.cut(4, 1));
    assert(!ett.connected(4, 0));
    assert(ett.link(4, 2));
    assert(ett.subtree_query(2, 0) == 27);          /// 10 + 17
    assert(ett.subtree_query(1, 0) == 24);          /// 13 + 11

    ett.set(0, 0);
    assert(ett.subtree_query(0, 2) == 24);

    return 0;
}
