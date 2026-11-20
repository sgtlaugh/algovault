#include "../common.h"

#define main library_main
#include "../../code_library/trees/leaf_peeling.cpp"
#undef main

/// BFS from root over adjacency lists: the reference parents, empty when the edges are not a tree on n nodes
vector<int> bfs_parents(int n, const vector<pair<int, int>>& edges, int root){
    if ((int)edges.size() != n - 1) return {};
    vector<vector<int>> adj(n);
    for (auto [u, v] : edges){
        if (u == v) return {};
        adj[u].push_back(v), adj[v].push_back(u);
    }

    vector<int> parent(n, -2), queue = {root};
    parent[root] = -1;
    for (int i = 0; i < (int)queue.size(); i++){
        for (int v : adj[queue[i]]){
            if (parent[v] != -2) continue;
            parent[v] = queue[i];
            queue.push_back(v);
        }
    }

    if ((int)queue.size() != n) return {};
    return parent;
}

/// order is a permutation ending at root with every node after its children, parents match BFS
void check(int n, const vector<pair<int, int>>& edges, int root){
    vector<int> expected = bfs_parents(n, edges, root);
    assert(!expected.empty());
    PeeledTree t = peel_leaves(n, edges, root);
    assert(t.parent == expected);

    assert((int)t.order.size() == n && t.order.back() == root);
    vector<int> pos(n, -1);
    for (int i = 0; i < n; i++){
        assert(0 <= t.order[i] && t.order[i] < n && pos[t.order[i]] == -1);
        pos[t.order[i]] = i;
    }
    for (int v = 0; v < n; v++){
        if (v != root) assert(pos[v] < pos[t.parent[v]]);
    }
}

/// Random labeled tree with a chosen shape: random attachment, path, star, caterpillar, broom
vector<pair<int, int>> random_tree(int n){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    int shape = stress::rand_int(0, 4), spine = max(1, n / 2);
    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i < spine ? i - 1 : stress::rand_int(0, spine - 1);
        else p = i < spine ? i - 1 : spine - 1;
        if (stress::rand_int(0, 1)) edges.push_back({label[i], label[p]});
        else edges.push_back({label[p], label[i]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Every multiset of n - 1 edges over [0, n): trees match BFS from every root, the rest abort
void check_all_edge_multisets(int n){
    vector<pair<int, int>> pairs;
    for (int u = 0; u < n; u++){
        for (int v = u; v < n; v++) pairs.push_back({u, v});
    }

    int m = pairs.size();
    vector<int> pick(max(0, n - 1), 0);
    while (true){
        vector<pair<int, int>> edges;
        for (int i = 0; i < n - 1; i++){
            auto [u, v] = pairs[pick[i]];
            if (i % 2) swap(u, v);
            edges.push_back({u, v});
        }

        if (!bfs_parents(n, edges, 0).empty()){
            for (int root = 0; root < n; root++) check(n, edges, root);
        }
        else{
            int root = stress::rand_int(0, n - 1);  /// one root per non-tree keeps the forked abort checks few
            stress::expect_assert("trees/leaf_peeling.cpp", [&]{ peel_leaves(n, edges, root); });
        }

        int i = n - 2;
        while (i >= 0 && pick[i] == m - 1) i--;
        if (i < 0) break;
        pick[i]++;
        for (int j = i + 1; j < n - 1; j++) pick[j] = pick[i];
    }
}

void check_rejects(){
    for (int n = 1; n <= 4; n++) check_all_edge_multisets(n);

    vector<pair<int, vector<pair<int, int>>>> bad = {
        {3, {{0, 3}, {1, 2}}},
        {3, {{-1, 0}, {1, 2}}},
        {3, {{0, 1}}},
        {3, {{0, 1}, {1, 2}, {2, 0}}},
        {6, {{0, 1}, {2, 3}, {3, 4}, {4, 5}, {5, 2}}},  /// a cycle away from the root, every root degree fine
        {6, {{0, 1}, {1, 2}, {3, 4}, {4, 5}, {5, 3}}},  /// the root's part is a path that peels fully
        {7, {{1, 2}, {0, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 3}}},  /// a split-off edge whose two ends peel each other
        {0, {}},
    };
    for (auto& [n, edges] : bad){
        stress::expect_assert("trees/leaf_peeling.cpp", [&]{ peel_leaves(n, edges, 0); });
    }

    /// a doubled edge plus a split-off edge: without the leaf-degree check the peel wraps around to all 4 nodes
    stress::expect_assert("trees/leaf_peeling.cpp", []{ peel_leaves(4, {{0, 3}, {0, 3}, {1, 2}}, 3); });
    stress::expect_assert("trees/leaf_peeling.cpp", []{ peel_leaves(3, {{0, 1}, {1, 2}}, 3); });
    stress::expect_assert("trees/leaf_peeling.cpp", []{ peel_leaves(3, {{0, 1}, {1, 2}}, -1); });
}

int main(){
    check(1, {}, 0);
    check(2, {{1, 0}}, 0);
    check(2, {{1, 0}}, 1);

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 500 ? it % 60 + 1 : stress::rand_int(1, 300);
        check(n, random_tree(n), stress::rand_int(0, n - 1));
    }

    for (long long it = 0; it < stress::scaled(2); it++){
        int n = stress::rand_int(100000, 200000);
        check(n, random_tree(n), stress::rand_int(0, n - 1));
    }

    check_rejects();
    return 0;
}
