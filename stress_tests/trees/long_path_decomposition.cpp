#include "../common.h"
#include <sys/wait.h>
#include <unistd.h>

#define main library_main
#include "../../code_library/trees/long_path_decomposition.cpp"
#undef main

struct NaiveTree{
    int n;
    vector<int> parent, level, tin, tout;
    vector<vector<int>> by_depth;

    NaiveTree(int n, int root, const vector<pair<int, int>>& edges) : n(n), parent(n, -1), level(n, 0), tin(n), tout(n){
        vector<vector<int>> adj(n);
        for (auto [u, v] : edges) adj[u].push_back(v), adj[v].push_back(u);

        vector<pair<int, int>> stack = {{root, 0}};
        int timer = 0;
        while (!stack.empty()){
            auto [u, i] = stack.back();
            if (i == 0) tin[u] = timer++;
            if (i == (int)adj[u].size()){
                tout[u] = timer - 1;
                stack.pop_back();
                continue;
            }
            stack.back().second++;
            int v = adj[u][i];
            if (v == parent[u] || v == root) continue;
            parent[v] = u, level[v] = level[u] + 1;
            stack.push_back({v, 0});
        }

        by_depth.assign(*max_element(level.begin(), level.end()) + 1, {});
        for (int v = 0; v < n; v++) by_depth[level[v]].push_back(tin[v]);
        for (auto& tins : by_depth) sort(tins.begin(), tins.end());
    }

    int kth(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        while (k--) v = parent[v];
        return v;
    }

    /// Nodes of v's subtree at depth level[v] + d are the nodes of that depth whose entry time falls in v's range
    int count_below(int v, int d) const{
        int target = level[v] + d;
        if (target >= (int)by_depth.size()) return 0;
        const auto& tins = by_depth[target];
        return upper_bound(tins.begin(), tins.end(), tout[v]) - lower_bound(tins.begin(), tins.end(), tin[v]);
    }

    /// O(n * depth): every node adds itself to each ancestor at its distance
    vector<vector<int>> all_counts() const{
        vector<vector<int>> cnt(n);
        for (int u = 0; u < n; u++){
            for (int a = u, d = 0; a != -1; a = parent[a], d++){
                if ((int)cnt[a].size() <= d) cnt[a].resize(d + 1, 0);
                cnt[a][d]++;
            }
        }
        return cnt;
    }
};

/// Random tree of a given shape with random labels
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
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else if (shape == 4) p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        else if (shape == 5) p = (i - 1) / 2;
        else p = stress::rand_int(max(0, i - 3), i - 1);
        edges.push_back({label[i], label[p]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Small trees: every depth array against the O(n * depth) counts, every (v, k) ancestor pair
void check_small(LongPathDecomposition& tree, const NaiveTree& naive, int n){
    auto expected = naive.all_counts();
    vector<bool> seen(n, false);
    int calls = 0;
    tree.depth_counts([&](int v, const int* cnt, int len){
        assert(!seen[v] && (naive.parent[v] == -1 || !seen[naive.parent[v]]));
        seen[v] = true, calls++;
        assert(len == (int)expected[v].size() && tree.height(v) == len - 1);
        assert(equal(cnt, cnt + len, expected[v].begin()));
    });
    assert(calls == n);

    for (int v = 0; v < n; v++){
        assert(tree.depth(v) == naive.level[v]);
        for (int k = -2; k <= naive.level[v] + 2; k++) assert(tree.kth_ancestor(v, k) == naive.kth(v, k));
    }
}

/// Large trees: sampled depths against the Euler tour counts, random ancestor queries
void check_large(LongPathDecomposition& tree, const NaiveTree& naive, int n, int queries){
    vector<bool> seen(n, false);
    int calls = 0;
    tree.depth_counts([&](int v, const int* cnt, int len){
        assert(!seen[v] && (naive.parent[v] == -1 || !seen[naive.parent[v]]));
        seen[v] = true, calls++;
        assert(naive.count_below(v, len - 1) > 0 && naive.count_below(v, len) == 0);
        assert(cnt[0] == 1 && cnt[len - 1] == naive.count_below(v, len - 1));
        for (int s = 0; s < 3; s++){
            int d = stress::rand_int(0, len - 1);
            assert(cnt[d] == naive.count_below(v, d));
        }
    });
    assert(calls == n);

    for (int q = 0; q < queries; q++){
        int v = stress::rand_int(0, n - 1);
        assert(tree.depth(v) == naive.level[v]);
        int k = q % 4 == 0 ? naive.level[v] - stress::rand_int(0, min(naive.level[v], 3)) : stress::rand_int(-1, naive.level[v] + 1);
        assert(tree.kth_ancestor(v, k) == naive.kth(v, k));
    }
}

void check(int n, int shape, bool small, int queries = 0){
    auto edges = make_tree(n, shape);
    LongPathDecomposition tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);

    for (int rebuild = 0; rebuild < 2; rebuild++){
        int root = stress::rand_int(0, n - 1);
        tree.build(root);
        NaiveTree naive(n, root, edges);
        if (small) check_small(tree, naive, n);
        else check_large(tree, naive, n, queries);
    }
}

/// build must abort on a non-tree, the alarm turns an endless traversal into a failure instead of a hang
void check_rejects(int n, const vector<pair<int, int>>& edges){
    pid_t pid = fork();
    assert(pid >= 0);
    if (pid == 0){
        alarm(5);
        assert(freopen("/dev/null", "w", stderr));
        LongPathDecomposition tree(n);
        for (auto [u, v] : edges) tree.add_edge(u, v);
        tree.build(0);
        _exit(0);
    }

    int status;
    assert(waitpid(pid, &status, 0) == pid);
    assert(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT);
}

int main(){
    check_rejects(4, {{0, 1}, {1, 2}, {2, 0}});          /// cycle plus an isolated node, n - 1 edges
    check_rejects(5, {{0, 1}, {1, 2}, {2, 3}, {3, 1}});  /// cycle below the root plus an isolated node
    check_rejects(3, {{0, 1}, {1, 2}, {2, 0}});          /// cycle, n edges
    check_rejects(3, {{0, 1}, {0, 1}, {1, 2}});          /// multi-edge
    check_rejects(3, {{0, 1}});                          /// disconnected
    check_rejects(2, {{0, 0}, {0, 1}});                  /// self-loop at the root, only the edge count sees it

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 400 ? it % 40 + 1 : stress::rand_int(1, 150);
        check(n, it % 7, true);
    }

    for (long long it = 0; it < stress::scaled(14); it++){
        check(stress::rand_int(20000, 60000), it % 7, false, 3000);
    }

    check(200000, 1, false, 100);  /// a path, deep enough to overflow an 8 MB stack with recursive traversals
    check(200000, 4, false, 100);

    return 0;
}
