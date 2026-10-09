#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/link_cut_tree.cpp"
#undef main

/// Naive forest: adjacency sets, one marked root per tree, a BFS per query
struct NaiveForest{
    int n;
    vector<set<int>> adj;
    vector<char> is_root;
    vector<long long> value;

    NaiveForest(int n, long long initial) : n(n), adj(n), is_root(n, 1), value(n, initial) {}

    /// BFS from u, returns parent links (-2 for unreached, -1 for u)
    vector<int> bfs(int u) const{
        vector<int> par(n, -2), queue = {u};
        par[u] = -1;
        for (int i = 0; i < (int)queue.size(); i++){
            for (int w : adj[queue[i]]){
                if (par[w] == -2) par[w] = queue[i], queue.push_back(w);
            }
        }
        return par;
    }

    int find_root(int u) const{
        vector<int> par = bfs(u);
        for (int x = 0; x < n; x++){
            if (par[x] != -2 && is_root[x]) return x;
        }
        assert(false);
    }

    bool connected(int u, int v) const{
        return bfs(u)[v] != -2;
    }

    void make_root(int u){
        is_root[find_root(u)] = 0;
        is_root[u] = 1;
    }

    bool link(int u, int v){
        if (connected(u, v)) return false;
        is_root[find_root(u)] = 0;
        adj[u].insert(v), adj[v].insert(u);
        return true;
    }

    bool cut(int u, int v){
        if (!adj[u].count(v)) return false;
        adj[u].erase(v), adj[v].erase(u);

        vector<int> par = bfs(u);
        bool u_has_root = false;
        for (int x = 0; x < n; x++) u_has_root |= par[x] != -2 && is_root[x];
        is_root[u_has_root ? v : u] = 1;
        return true;
    }

    vector<int> path(int u, int v) const{
        vector<int> par = bfs(u), res;
        for (int x = v; x != -1; x = par[x]) res.push_back(x);
        return res;
    }

    int parent(int u) const{
        return bfs(find_root(u))[u];
    }

    int lca(int u, int v) const{
        if (!connected(u, v)) return -1;
        vector<int> par = bfs(find_root(u)), on_path(n, 0);
        for (int x = u; x != -1; x = par[x]) on_path[x] = 1;
        int x = v;
        while (!on_path[x]) x = par[x];
        return x;
    }
};

template <typename Combine>
void check_random(int n, int ops, long long identity, Combine combine){
    LinkCutTree<long long, Combine> lct(n, identity, combine);
    NaiveForest ref(n, identity);
    int max_value = stress::rand_int(0, 1) ? 3 : 1000000000;

    for (int op = 0; op < ops; op++){
        int kind = stress::rand_int(0, 10), u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (kind <= 2) assert(lct.link(u, v) == ref.link(u, v));
        else if (kind == 3){
            if (!ref.adj[u].empty() && stress::rand_int(0, 3)) v = *ref.adj[u].begin();
            assert(lct.cut(u, v) == ref.cut(u, v));
        }
        else if (kind == 4) lct.make_root(u), ref.make_root(u);
        else if (kind == 5){
            long long x = stress::rand_int(-max_value, max_value);
            lct.set(u, x), ref.value[u] = x;
        }
        else if (kind == 6) assert(lct.connected(u, v) == ref.connected(u, v));
        else if (kind == 7) assert(lct.find_root(u) == ref.find_root(u));
        else if (kind == 8) assert(lct.lca(u, v) == ref.lca(u, v));
        else if (kind == 9) assert(lct.parent(u) == ref.parent(u));
        else if (ref.connected(u, v)){
            long long expected = identity;
            for (int x : ref.path(u, v)) expected = combine(expected, ref.value[x]);
            assert(lct.query(u, v) == expected);
            assert(lct.find_root(u) == ref.find_root(u));
        }
    }

    for (int u = 0; u < n; u++) assert(lct.find_root(u) == ref.find_root(u));
}

/// A path 0-1-...-n-1 linked child to parent: maximal depth, checked against prefix sums
void check_long_path(int n){
    LinkCutTree<long long> lct(n);
    vector<long long> prefix(n + 1, 0);
    for (int i = 0; i < n; i++){
        long long x = stress::rand_int(-1000000000, 1000000000);
        lct.set(i, x);
        prefix[i + 1] = prefix[i] + x;
    }

    for (int i = 1; i < n; i++) assert(lct.link(i, i - 1));
    for (int it = 0; it < 20000; it++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(0, n - 1);
        assert(lct.lca(l, r) == min(l, r));
        assert(lct.query(l, r) == prefix[max(l, r) + 1] - prefix[min(l, r)]);
    }

    /// Sequential access patterns expose a broken splay rule (move-to-root) as a timeout, random ones do not
    for (int r = 0; r < 3; r++){
        for (int i = 0; i < n; i++) assert(lct.lca(i, n - 1 - i) == min(i, n - 1 - i));
    }
    for (int i = 0; i < n; i++) assert(lct.parent(i) == i - 1);

    int mid = n / 2;
    assert(lct.cut(mid, mid - 1) && !lct.connected(0, n - 1));
    assert(lct.find_root(n - 1) == mid && lct.find_root(mid - 1) == 0);
    lct.make_root(n - 1);
    assert(lct.lca(mid, n - 2) == n - 2 && lct.query(n - 1, mid) == prefix[n] - prefix[mid]);
}

int main(){
    auto add = [](long long a, long long b){ return a + b; };
    auto mx = [](long long a, long long b){ return max(a, b); };

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it % 8 + 1 : stress::rand_int(1, 60);
        int ops = stress::rand_int(1, 300);
        if (it % 2) check_random(n, ops, 0LL, add);
        else check_random(n, ops, LLONG_MIN, mx);
    }

    check_long_path(200000);

    return 0;
}
