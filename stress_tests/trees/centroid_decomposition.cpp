#include "../common.h"

#define main library_main
#include "../../code_library/trees/centroid_decomposition.cpp"
#undef main

vector<tuple<int, int, long long>> make_forest(int n, int shape, long long max_w, bool forest){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<tuple<int, int, long long>> edges;
    for (int i = 1; i < n; i++){
        if (forest && stress::rand_int(0, 4) == 0) continue;
        int p = shape == 0 ? stress::rand_int(0, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 0 : (i - 1) / 2;
        edges.push_back({label[i], label[p], stress::rand_int(1, max_w)});
    }

    return edges;
}

/// All-pairs distances by BFS from every node, -1 between different trees
vector<vector<long long>> all_distances(int n, const vector<tuple<int, int, long long>>& edges){
    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges) adj[u].push_back({v, w}), adj[v].push_back({u, w});

    vector<vector<long long>> dist(n, vector<long long>(n, -1));
    for (int s = 0; s < n; s++){
        vector<int> queue = {s};
        dist[s][s] = 0;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (auto [v, w] : adj[u]){
                if (dist[s][v] == -1) dist[s][v] = dist[s][u] + w, queue.push_back(v);
            }
        }
    }

    return dist;
}

/// Path counts against all-pairs BFS; branches against the centroid tree and true distances; centroid sizes halve
void check(int n, int shape, long long max_w, bool forest){
    auto edges = make_forest(n, shape, max_w, forest);
    auto dist = all_distances(n, edges);
    CentroidDecomposition cd(n);
    for (auto [u, v, w] : edges) cd.add_edge(u, v, w);

    vector<CentroidDecomposition::Branches> seen(n);
    vector<int> visit_order, done(n, 0);
    cd.build([&](int c, const CentroidDecomposition::Branches& branches){
        assert(!done[c]);
        done[c] = 1, seen[c] = branches, visit_order.push_back(c);
    });
    assert((int)visit_order.size() == n);

    vector<vector<int>> children(n);
    int roots = 0;
    for (int c : visit_order){
        if (cd.parent[c] == -1){
            roots++;
            assert(cd.depth[c] == 0);
            continue;
        }
        assert(done[cd.parent[c]] && cd.depth[c] == cd.depth[cd.parent[c]] + 1);
        children[cd.parent[c]].push_back(c);
    }
    assert(roots == n - (int)edges.size());

    vector<int> pos(n);
    for (int i = 0; i < n; i++) pos[visit_order[i]] = i;
    vector<vector<int>> subtree(n);
    for (int i = n - 1; i >= 0; i--){
        int c = visit_order[i];
        subtree[c].push_back(c);
        for (int d : children[c]){
            assert(pos[d] > pos[c]);
            subtree[c].insert(subtree[c].end(), subtree[d].begin(), subtree[d].end());
        }
        sort(subtree[c].begin(), subtree[c].end());
    }

    for (int c = 0; c < n; c++){
        assert(cd.depth[c] <= __lg(n));
        vector<vector<int>> pieces;
        for (auto& branch : seen[c]){
            vector<int> nodes;
            for (auto [v, d] : branch){
                assert(dist[c][v] == d);
                nodes.push_back(v);
            }
            sort(nodes.begin(), nodes.end());
            pieces.push_back(nodes);
        }
        vector<vector<int>> expected;
        for (int d : children[c]){
            assert(2 * subtree[d].size() <= subtree[c].size());
            expected.push_back(subtree[d]);
        }
        sort(pieces.begin(), pieces.end()), sort(expected.begin(), expected.end());
        assert(pieces == expected);
    }

    for (long long k : {0LL, 1LL, (long long)stress::rand_int(0, 3LL * max_w), (long long)stress::rand_int(0, (long long)n * max_w), (long long)n * max_w}){
        long long brute = 0;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) brute += dist[u][v] != -1 && dist[u][v] <= k;
        }
        assert(count_paths_at_most(cd, k) == brute);
    }

    auto parent = cd.parent, depth = cd.depth;
    cd.build();
    assert(cd.parent == parent && cd.depth == depth);
}

/// build must abort on a non-forest instead of looping forever around a cycle
void check_rejects(int n, const vector<pair<int, int>>& edges){
    stress::expect_assert("trees/centroid_decomposition.cpp", [&]{
        CentroidDecomposition cd(n);
        for (auto [u, v] : edges) cd.add_edge(u, v);
        cd.build();
    });
}

int main(){
    check_rejects(4, {{0, 1}, {1, 2}, {2, 0}});                  /// cycle plus an isolated node, n - 1 edges
    check_rejects(6, {{0, 1}, {2, 3}, {3, 4}, {4, 5}, {5, 3}});  /// cycle in the second tree, below its first node
    check_rejects(3, {{0, 1}, {0, 1}, {1, 2}});                  /// multi-edge
    check_rejects(2, {{0, 0}, {0, 1}});                          /// self-loop at the first node
    check_rejects(4, {{0, 1}, {1, 2}, {2, 2}, {2, 3}});          /// self-loop in the middle of a path

    for (long long it = 0; it < stress::scaled(700); it++){
        int n = it < 200 ? it % 12 + 1 : stress::rand_int(1, 40);
        check(n, it % 4, it % 3 ? 1 : 20, it % 5 == 0);
    }
    check(1500, 0, 1, false);
    check(1500, 1, 1000000000, false);
    check(1500, 0, 4000000000000LL, false);  /// weights beyond int, so a weight narrowed anywhere breaks the distances
    check(1024, 1, 1, false);
    check(1500, 2, 5, true);

    /// A path of 100000 nodes: deep enough to overflow an 8 MB stack with recursion, and a pair count beyond int
    int n = 100000;
    CentroidDecomposition path(n);
    for (int i = 1; i < n; i++) path.add_edge(i - 1, i);
    long long delivered = 0;
    path.build([&](int, const CentroidDecomposition::Branches& branches){
        for (auto& branch : branches) delivered += branch.size();
    });
    assert(*max_element(path.depth.begin(), path.depth.end()) <= __lg(n));
    assert(delivered <= (long long)n * __lg(n));
    assert(count_paths_at_most(path, n) == (long long)n * (n - 1) / 2);

    return 0;
}
