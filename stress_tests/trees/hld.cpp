#include "../common.h"

#define main library_main
#include "../../code_library/trees/hld.cpp"
#undef main

vector<pair<int, int>> make_tree(int n, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = (i - 1) / 2;
        else p = i < n / 2 ? i - 1 : stress::rand_int(0, i - 1);
        edges.push_back({label[i], label[p]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Path ranges, updates, subtrees and lca against parent climbing on the same rooted tree
void check(int n, int shape, int ops){
    auto edges = make_tree(n, shape);
    int root = stress::rand_int(0, n - 1);
    HLD hld(n);
    for (auto [u, v] : edges) hld.add_edge(u, v);
    hld.build(root);

    /// Own BFS from the chosen root, so a build that roots the tree elsewhere cannot agree with itself
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges) adj[u].push_back(v), adj[v].push_back(u);
    vector<int> parent(n, -1), depth(n, -1), bfs = {root};
    depth[root] = 0;
    for (int i = 0; i < (int)bfs.size(); i++){
        for (int v : adj[bfs[i]]){
            if (depth[v] == -1) parent[v] = bfs[i], depth[v] = depth[bfs[i]] + 1, bfs.push_back(v);
        }
    }
    assert(hld.parent == parent && hld.depth == depth);

    vector<int> seen(n, 0);
    for (int v = 0; v < n; v++) seen[hld.pos[v]]++;
    for (int i = 0; i < n; i++) assert(seen[i] == 1);

    auto naive_path = [&](int u, int v){
        vector<int> left, right;
        while (u != v){
            if (depth[u] >= depth[v]) left.push_back(u), u = parent[u];
            else right.push_back(v), v = parent[v];
        }
        left.push_back(u);
        left.insert(left.end(), right.rbegin(), right.rend());
        return make_pair(left, u);
    };

    vector<long long> value(n), at(n);
    for (int v = 0; v < n; v++) value[v] = stress::rand_int(-1000000000, 1000000000), at[hld.pos[v]] = value[v];

    for (int op = 0; op < ops; op++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        bool edge_mode = stress::rand_int(0, 1);
        auto [nodes, l] = naive_path(u, v);
        assert(hld.lca(u, v) == l);

        vector<pair<int, int>> ranges = hld.path(u, v, edge_mode);
        vector<int> covered;
        for (auto [a, b] : ranges){
            assert(a <= b);
            for (int i = a; i <= b; i++) covered.push_back(i);
        }
        vector<int> expected;
        for (int x : nodes){
            if (!(edge_mode && x == l)) expected.push_back(hld.pos[x]);
        }
        sort(covered.begin(), covered.end());
        sort(expected.begin(), expected.end());
        assert(covered == expected);

        if (op % 3 == 0){
            long long delta = stress::rand_int(-1000, 1000);
            for (auto [a, b] : ranges){
                for (int i = a; i <= b; i++) at[i] += delta;
            }
            for (int x : nodes){
                if (!(edge_mode && x == l)) value[x] += delta;
            }
        }
        long long fast = 0, slow = 0;
        for (auto [a, b] : hld.path(u, v, edge_mode)){
            for (int i = a; i <= b; i++) fast += at[i];
        }
        for (int x : nodes){
            if (!(edge_mode && x == l)) slow += value[x];
        }
        assert(fast == slow);

        auto [a, b] = hld.subtree(u);
        vector<int> in_range, in_subtree;
        for (int i = a; i <= b; i++) in_range.push_back(i);
        for (int x = 0; x < n; x++){
            int y = x;
            while (y != -1 && y != u) y = parent[y];
            if (y == u) in_subtree.push_back(hld.pos[x]);
        }
        sort(in_subtree.begin(), in_subtree.end());
        if (op % 7 == 0) assert(in_range == in_subtree);
    }
}

int main(){
    /// Rebuilding from another root must match a fresh decomposition
    for (long long it = 0; it < stress::scaled(500); it++){
        int n = stress::rand_int(1, 40);
        auto edges = make_tree(n, it % 5);
        HLD reused(n);
        for (auto [u, v] : edges) reused.add_edge(u, v);
        for (int round = 0; round < 3; round++){
            int root = stress::rand_int(0, n - 1);
            HLD fresh(n);
            for (auto [u, v] : edges) fresh.add_edge(u, v);
            reused.build(root), fresh.build(root);
            assert(reused.pos == fresh.pos && reused.head == fresh.head && reused.size == fresh.size);
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 200 ? it % 30 + 1 : stress::rand_int(1, 120);
        check(n, it % 5, 60);
    }

    for (long long it = 0; it < stress::scaled(5); it++){
        int n = 200000;
        auto edges = make_tree(n, it % 5);
        HLD hld(n);
        for (auto [u, v] : edges) hld.add_edge(u, v);
        hld.build(stress::rand_int(0, n - 1));
        for (int q = 0; q < 2000; q++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            assert(hld.path(u, v).size() <= 2 * 18 + 2);
        }
    }

    return 0;
}
