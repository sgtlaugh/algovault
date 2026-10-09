#include "../common.h"

#define main library_main
#include "../../code_library/graphs/three_edge_connected_components.cpp"
#undef main

/// Component label of every node using all edges except edges a and b
vector<int> labels(int n, const vector<pair<int, int>>& edges, int a, int b){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };

    for (int i = 0; i < (int)edges.size(); i++){
        if (i != a && i != b) parent[find(edges[i].first)] = find(edges[i].second);
    }

    vector<int> res(n);
    for (int i = 0; i < n; i++) res[i] = find(i);
    return res;
}

/// u and v share a class when every removal of at most 2 edges leaves them connected, groups ordered by smallest vertex
vector<vector<int>> brute(int n, const vector<pair<int, int>>& edges){
    int m = edges.size();
    vector<vector<char>> together(n, vector<char>(n, 1));
    for (int a = -1; a < m; a++){
        for (int b = a; b < m; b++){
            auto lab = labels(n, edges, a, b);
            for (int u = 0; u < n; u++){
                for (int v = u + 1; v < n; v++){
                    if (lab[u] != lab[v]) together[u][v] = 0;
                }
            }
        }
    }

    vector<vector<int>> groups;
    vector<char> used(n, 0);
    for (int u = 0; u < n; u++){
        if (used[u]) continue;
        groups.push_back({});
        for (int v = u; v < n; v++){
            if (together[u][v]) used[v] = 1, groups.back().push_back(v);
        }
    }
    return groups;
}

/// Random multigraphs: dense, sparse, or cycles glued together, with parallel edges and self loops
vector<pair<int, int>> random_graph(int n){
    vector<pair<int, int>> edges;
    int shape = stress::rand_int(0, 2);

    if (shape == 0){
        int density = stress::rand_int(0, n > 9 ? 35 : 100);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (stress::rand_int(0, 99) < density) edges.push_back({u, v});
            }
        }
    }
    else if (shape == 1){
        for (int c = stress::rand_int(0, 2 * n); c; c--) edges.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1)});
    }
    else{
        for (int c = stress::rand_int(1, n); c; c--){
            vector<int> cycle(stress::rand_int(1, min(n, 5)));
            for (auto& x : cycle) x = stress::rand_int(0, n - 1);
            for (size_t i = 0; i < cycle.size(); i++) edges.push_back({cycle[i], cycle[(i + 1) % cycle.size()]});
        }
    }

    for (int extra = stress::rand_int(0, 3); extra && !edges.empty(); extra--) edges.push_back(edges[stress::rand_int(0, edges.size() - 1)]);
    if (stress::rand_int(0, 1)){
        int x = stress::rand_int(0, n - 1);
        edges.push_back({x, x});
    }
    for (auto& [u, v] : edges){
        if (stress::rand_int(0, 1)) swap(u, v);
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

/// Checks the groups and the comp labels against each other
void check(ThreeEdgeConnected& g, const vector<vector<int>>& groups, const vector<vector<int>>& expected){
    assert(groups == expected);
    assert(g.count == (int)expected.size());
    for (int c = 0; c < g.count; c++){
        for (int v : expected[c]) assert(g.comp[v] == c);
    }
}

/// k blocks of K4 in a chain or ring, consecutive blocks joined by `links` edges between random endpoints
/// A chain with links <= 2 splits into the blocks, links >= 3 or a ring with links >= 2 merges everything
void blocks_test(int k, int links, bool ring){
    int n = 4 * k;
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    ThreeEdgeConnected g(n);
    for (int b = 0; b < k; b++){
        for (int i = 0; i < 4; i++){
            for (int j = i + 1; j < 4; j++) g.add_edge(label[4 * b + i], label[4 * b + j]);
        }
    }
    for (int b = 0; b + 1 < k || (ring && k > 2 && b < k); b++){
        int c = (b + 1) % k;
        for (int l = 0; l < links; l++) g.add_edge(label[4 * b + stress::rand_int(0, 3)], label[4 * c + stress::rand_int(0, 3)]);
    }

    auto start = chrono::steady_clock::now();
    auto groups = g.run();
    assert(chrono::steady_clock::now() - start < chrono::seconds(3));

    bool merged = links >= 3 || (ring && k > 2 && links >= 2);
    assert(g.count == (merged ? 1 : k));
    for (int b = 0; b < k; b++){
        for (int i = 0; i < 4; i++) assert(g.comp[label[4 * b + i]] == (merged ? 0 : g.comp[label[4 * b]]));
    }
    for (auto& group : groups) assert(is_sorted(group.begin(), group.end()));
}

/// Path over shuffled labels where each step carries `copies` parallel edges, so m = copies * (n - 1) dominates n
/// Cutting one step needs all of its copies removed: copies <= 2 isolates every vertex, copies >= 3 keeps one class
void parallel_path_test(int n, int copies){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<pair<int, int>> edges;
    for (int i = 0; i + 1 < n; i++){
        for (int c = 0; c < copies; c++) edges.push_back({label[i], label[i + 1]});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    ThreeEdgeConnected g(n);
    for (auto [u, v] : edges) g.add_edge(u, v);

    auto start = chrono::steady_clock::now();
    g.run();
    assert(chrono::steady_clock::now() - start < chrono::seconds(3));

    assert(g.count == (copies >= 3 ? 1 : n));
}

int main(){
    /// Every multigraph on 4 vertices with up to 5 edges drawn from the 10 vertex pairs including loops
    vector<pair<int, int>> pairs;
    for (int u = 0; u < 4; u++){
        for (int v = u; v < 4; v++) pairs.push_back({u, v});
    }
    for (int m = 0; m <= 5; m++){
        vector<int> pick(m, 0);
        while (true){
            vector<pair<int, int>> edges;
            for (int p : pick) edges.push_back(pairs[p]);
            ThreeEdgeConnected g(4);
            for (auto [u, v] : edges) g.add_edge(u, v);
            auto groups = g.run();
            check(g, groups, brute(4, edges));

            int i = m - 1;
            while (i >= 0 && pick[i] == (int)pairs.size() - 1) i--;
            if (i < 0) break;
            pick[i]++;
            for (int j = i + 1; j < m; j++) pick[j] = pick[i];
        }
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, it % 20 ? 9 : 16);
        auto edges = random_graph(n);
        ThreeEdgeConnected g(n);
        for (auto [u, v] : edges) g.add_edge(u, v);

        auto groups = g.run();
        check(g, groups, brute(n, edges));

        /// run() again after more edges must not read state left by the first run
        if (it % 7 == 0){
            edges.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1)});
            g.add_edge(edges.back().first, edges.back().second);
            groups = g.run();
            check(g, groups, brute(n, edges));
        }
    }

    /// Deep DFS: 25000 blocks in a chain is a 100000 vertex graph, deeper than a recursive DFS survives under ASan with 8 MB
    for (int links : {1, 2, 3}){
        blocks_test(25000, links, false);
        blocks_test(25000, links, true);
    }
    for (int k : {1, 2, 3}){
        for (int links : {0, 1, 2, 3}){
            blocks_test(k, links, false);
            blocks_test(k, links, true);
        }
    }

    parallel_path_test(200000, 2);
    parallel_path_test(200000, 5);

    ThreeEdgeConnected path(300000);
    for (int i = 0; i + 1 < 300000; i++) path.add_edge(i, i + 1);
    path.run();
    assert(path.count == 300000);
    return 0;
}
