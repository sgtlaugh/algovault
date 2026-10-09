#include "../common.h"

#define main library_main
#include "../../code_library/trees/virtual_tree.cpp"
#undef main

struct NaiveTree{
    vector<int> parent, level, order;
    vector<long long> up_weight;

    NaiveTree(int n, int root, const vector<array<long long, 3>>& edges) : parent(n, -1), level(n, 0), up_weight(n, 0){
        vector<vector<pair<int, long long>>> adj(n);
        for (auto [u, v, w] : edges) adj[u].push_back({(int)v, w}), adj[v].push_back({(int)u, w});

        vector<bool> seen(n, false);
        order = {root}, seen[root] = true;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (auto [v, w] : adj[u]){
                if (seen[v]) continue;
                seen[v] = true, parent[v] = u, level[v] = level[u] + 1, up_weight[v] = w;
                order.push_back(v);
            }
        }
    }

    int lca(int u, int v) const{
        while (level[u] > level[v]) u = parent[u];
        while (level[v] > level[u]) v = parent[v];
        while (u != v) u = parent[u], v = parent[v];
        return u;
    }

    long long dist(int u, int v) const{
        long long res = 0;
        while (level[u] > level[v]) res += up_weight[u], u = parent[u];
        while (level[v] > level[u]) res += up_weight[v], v = parent[v];
        while (u != v) res += up_weight[u] + up_weight[v], u = parent[u], v = parent[v];
        return res;
    }
};

vector<array<long long, 3>> make_tree(int n, int shape, long long max_w){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<array<long long, 3>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        edges.push_back({label[i], label[p], stress::rand_int(0, max_w)});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Shared checks given the expected vertex set (as a mask) and each member's expected virtual parent
void check_result(const vector<pair<int, int>>& res, const vector<char>& in_set, const vector<int>& expected_parent, int distinct){
    int members = count(in_set.begin(), in_set.end(), 1);
    assert((int)res.size() == members);
    assert(distinct == 0 ? res.empty() : (int)res.size() <= 2 * distinct - 1);
    if (res.empty()) return;

    vector<int> stack;
    vector<char> seen(in_set.size(), 0);
    for (int i = 0; i < (int)res.size(); i++){
        auto [v, p] = res[i];
        assert(in_set[v] && !seen[v]);
        seen[v] = 1;
        assert(i == 0 ? p == -1 : (p >= 0 && p < i));
        assert(expected_parent[v] == (p == -1 ? -1 : res[p].first));

        /// DFS order: the parent must be on the current root path
        while (!stack.empty() && stack.back() != p) stack.pop_back();
        assert(i == 0 || !stack.empty());
        stack.push_back(i);
    }
}

/// Reference: close the set under pairwise naive lca until stable, parent = nearest proper ancestor in the closure
void check_small(const VirtualTree& tree, const NaiveTree& naive, const vector<int>& nodes){
    int n = naive.parent.size();
    vector<char> in_set(n, 0);
    vector<int> members;
    for (int v : nodes){
        if (!in_set[v]) in_set[v] = 1, members.push_back(v);
    }
    int distinct = members.size();

    for (bool changed = true; changed; ){
        changed = false;
        for (int i = 0; i < (int)members.size(); i++){
            for (int j = i + 1; j < (int)members.size(); j++){
                int l = naive.lca(members[i], members[j]);
                if (!in_set[l]) in_set[l] = 1, members.push_back(l), changed = true;
            }
        }
    }

    vector<int> expected_parent(n, -2);
    for (int v : members){
        int u = naive.parent[v];
        while (u != -1 && !in_set[u]) u = naive.parent[u];
        expected_parent[v] = u;
    }

    auto res = tree.compress(nodes);
    check_result(res, in_set, expected_parent, distinct);
    for (int i = 0; i < (int)res.size(); i++){
        auto [v, p] = res[i];
        assert(tree.depth(v) == naive.level[v]);
        if (p != -1) assert(tree.dist(v, res[p].first) == naive.dist(v, res[p].first));
    }
}

/// Reference in O(n): a vertex is kept iff it is chosen or at least two of its child subtrees hold chosen vertices
void check_large(const VirtualTree& tree, const NaiveTree& naive, const vector<int>& nodes){
    int n = naive.parent.size();
    vector<char> chosen(n, 0), in_set(n, 0);
    for (int v : nodes) chosen[v] = 1;
    int distinct = count(chosen.begin(), chosen.end(), 1);

    vector<int> below(n, 0), busy_children(n, 0);
    for (int i = n - 1; i >= 0; i--){
        int v = naive.order[i];
        below[v] += chosen[v];
        in_set[v] = chosen[v] || busy_children[v] >= 2;
        int p = naive.parent[v];
        if (p != -1 && below[v] > 0) below[p] += below[v], busy_children[p]++;
    }

    vector<int> nearest(n, -1), expected_parent(n, -2);
    for (int v : naive.order){
        int p = naive.parent[v];
        if (p != -1) nearest[v] = in_set[p] ? p : nearest[p];
        if (in_set[v]) expected_parent[v] = nearest[v];
    }

    check_result(tree.compress(nodes), in_set, expected_parent, distinct);
}

vector<int> random_nodes(int n, int k){
    vector<int> nodes(k);
    for (auto& v : nodes) v = stress::rand_int(0, n - 1);
    return nodes;
}

int main(){
    for (int n = 1; n <= 7; n++){
        for (int it = 0; it < 15; it++){
            auto edges = make_tree(n, it % 5, 9);
            int root = stress::rand_int(0, n - 1);
            VirtualTree tree(n);
            for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
            tree.build(root);
            NaiveTree naive(n, root, edges);

            for (int mask = 0; mask < (1 << n); mask++){
                vector<int> nodes;
                for (int v = 0; v < n; v++){
                    if (mask >> v & 1) nodes.push_back(v);
                }
                shuffle(nodes.begin(), nodes.end(), stress::rng());
                check_small(tree, naive, nodes);
                check_large(tree, naive, nodes);
            }
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = it < 300 ? it % 60 + 1 : stress::rand_int(1, 400);
        auto edges = make_tree(n, it % 5, it % 3 ? 1000000000LL : 1);
        VirtualTree tree(n);
        for (auto [u, v, w] : edges) tree.add_edge(u, v, w);

        /// the second pass rebuilds the same object at a new root, so stale data from the first build would surface
        for (int pass = 0; pass < 2; pass++){
            int root = stress::rand_int(0, n - 1);
            tree.build(root);
            NaiveTree naive(n, root, edges);

            for (int q = 0; q < 10; q++){
                auto nodes = random_nodes(n, stress::rand_int(0, min(n + 5, 40)));
                check_small(tree, naive, nodes);
                check_large(tree, naive, nodes);
            }
        }
    }

    for (long long it = 0; it < stress::scaled(10); it++){
        int n = stress::rand_int(20000, 60000);
        auto edges = make_tree(n, it % 5, 1000000000LL);
        int root = stress::rand_int(0, n - 1);
        VirtualTree tree(n);
        for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
        tree.build(root);
        NaiveTree naive(n, root, edges);

        for (int k : {1, 2, 10, 1000, n / 2, n, 2 * n}) check_large(tree, naive, random_nodes(n, k));
    }

    /// a path, deep enough to overflow an 8 MB stack with recursive traversals
    int n = 200000;
    auto edges = make_tree(n, 1, 1000000000LL);
    VirtualTree tree(n);
    for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
    tree.build(0);
    NaiveTree naive(n, 0, edges);
    for (int k : {2, 100, n}) check_large(tree, naive, random_nodes(n, k));

    return 0;
}
