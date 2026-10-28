#include "../common.h"

#define main library_main
#include "../../code_library/graphs/online_dynamic_connectivity.cpp"
#undef main

/// BFS component labels over the live edge copies
vector<int> bfs_labels(int n, const vector<pair<int, int>>& live){
    vector<vector<int>> adj(n);
    for (auto [u, v]: live){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> label(n, -1);
    for (int s = 0; s < n; s++){
        if (label[s] != -1) continue;
        queue<int> bfs;
        bfs.push(s), label[s] = s;
        while (!bfs.empty()){
            int u = bfs.front();
            bfs.pop();
            for (int v: adj[u]){
                if (label[v] == -1) label[v] = s, bfs.push(v);
            }
        }
    }

    return label;
}

int max_level_seen = 0;

void check(OnlineDynamicConnectivity& g, int n, const vector<pair<int, int>>& live, const vector<long long>& value){
    auto label = bfs_labels(n, live);
    vector<int> size(n, 0);
    vector<long long> sum(n, 0);
    int roots = 0;
    for (int v = 0; v < n; v++) size[label[v]]++, sum[label[v]] += value[v], roots += label[v] == v;

    assert(g.components == roots);
    for (int v = 0; v < n; v++){
        assert(g.size(v) == size[label[v]]);
        assert(g.sum(v) == sum[label[v]]);
    }
    for (int k = 0; k < min(n * n, 30); k++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        assert(g.connected(u, v) == (label[u] == label[v]));
    }

    /// Level i trees have at most n / 2^i vertices, the bound behind the O(log^2 n) amortized updates
    for (auto& e: g.edges){
        assert(e.level <= __lg(max(n, 1)));
        max_level_seen = max(max_level_seen, e.level);
    }
    for (int i = 0; i < (int)g.vertex_at.size(); i++){
        for (int x: g.vertex_at[i]){
            if (x == -1) continue;
            g.splay(x);
            assert(g.t[x].size <= n >> i);
        }
    }

    /// Catches a non-tree edge raised above the level that joins its endpoints, which answers hide for a while
    for (auto& e: g.edges){
        if (e.level == -1 || e.u == e.v) continue;
        if (!e.arcs.empty()) assert((int)e.arcs.size() == e.level + 1);
        else assert(g.same_tree(g.vertex_at[e.level][e.u], g.vertex_at[e.level][e.v]));
    }
}

/// Random operations, checked against BFS after every one when check_every is 1
void run(int n, int ops, int check_every){
    vector<long long> value(n);
    for (auto& x: value) x = stress::rand_int(-1000000000, 1000000000);
    OnlineDynamicConnectivity g(value);
    vector<pair<int, int>> live, removed;

    for (int i = 0; i < ops; i++){
        int type = stress::rand_int(0, 9);
        if (n == 0) break;

        if (type <= 3){
            /// Mostly short hops so cycles form, then removals need replacements; also re-adds and self loops
            int u = stress::rand_int(0, n - 1);
            int v = stress::rand_int(0, 2) ? (u + stress::rand_int(1, 3)) % n : stress::rand_int(0, n - 1);
            if (!removed.empty() && stress::rand_int(0, 3) == 0){
                tie(u, v) = removed[stress::rand_int(0, removed.size() - 1)];
            }
            g.add_edge(u, v);
            live.push_back({u, v});
        }
        else if (type <= 7 && !live.empty()){
            int k = stress::rand_int(0, live.size() - 1);
            auto [u, v] = live[k];
            swap(live[k], live.back());
            live.pop_back();
            removed.push_back({u, v});
            if (stress::rand_int(0, 1)) swap(u, v);
            g.remove_edge(u, v);
        }
        else{
            int v = stress::rand_int(0, n - 1);
            long long x = stress::rand_int(-1000000000, 1000000000);
            g.add(v, x);
            value[v] += x;
        }

        if (i % check_every == 0) check(g, n, live, value);
    }
    check(g, n, live, value);
}

int main(){
    for (int n = 0; n <= 4; n++) run(n, 200, 1);
    for (long long it = 0; it < stress::scaled(700); it++) run(stress::rand_int(1, it % 10 ? 12 : 40), 250, 1);
    assert(max_level_seen >= 2);  /// the replacement search above level 0 really ran, not only its level 0 pass

    for (int n: {300, 3000}){
        clock_t start = clock();
        run(n, 60000, n / 10);
        assert(clock() - start < 20 * CLOCKS_PER_SEC);
    }

    /// Two dense halves joined by a bridge that is deleted and re-added: each deletion scans a whole half for nothing
    int half = 400, rounds = 3000;
    OnlineDynamicConnectivity g(2 * half);
    for (int side = 0; side < 2; side++){
        for (int i = 0; i < half; i++){
            for (int d = 1; d <= 20; d++) g.add_edge(side * half + i, side * half + (i + d) % half);
        }
    }
    clock_t start = clock();
    for (int r = 0; r < rounds; r++){
        int u = stress::rand_int(0, half - 1), v = half + stress::rand_int(0, half - 1);
        g.add_edge(u, v);
        assert(g.components == 1);
        g.remove_edge(v, u);
        assert(g.components == 2);
        assert(g.size(u) == half);
    }
    assert(clock() - start < 10 * CLOCKS_PER_SEC);

    /// A path cut and relinked at random points: deep splay trees, every cut a tree edge with no replacement
    int len = 30000;
    OnlineDynamicConnectivity path(len);
    for (int i = 0; i + 1 < len; i++) path.add_edge(i, i + 1);
    assert(path.size(0) == len);
    for (int r = 0; r < 6000; r++){
        int i = stress::rand_int(0, len - 2);
        path.remove_edge(i, i + 1);
        assert(path.size(0) == i + 1);
        assert(!path.connected(i, i + 1));
        path.add_edge(i + 1, i);
    }
    assert(path.components == 1);

    /// Vertices past 46341, where an int edge key min * n + max overflows
    int big = 100000;
    OnlineDynamicConnectivity wide(big);
    wide.add_edge(big - 1, big - 2);
    wide.add_edge(big - 3, big - 1);
    wide.remove_edge(big - 2, big - 1);
    assert(wide.connected(big - 3, big - 1));
    assert(!wide.connected(big - 2, big - 1));
    assert(wide.components == big - 1);

    return 0;
}
