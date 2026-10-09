#include "../common.h"

#define main library_main
#include "../../code_library/graphs/online_bridges.cpp"
#undef main

/// Component label of every vertex using all edges except the skipped edge indices
vector<int> labels(int n, const vector<pair<int, int>>& edges, const vector<char>& skipped){
    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };

    for (int i = 0; i < (int)edges.size(); i++){
        if (!skipped[i]) parent[find(edges[i].first)] = find(edges[i].second);
    }

    vector<int> res(n);
    for (int i = 0; i < n; i++) res[i] = find(i);
    return res;
}

/// Bridge flags by removing each edge in turn, O(m^2)
vector<char> brute_bridges(int n, const vector<pair<int, int>>& edges){
    int m = edges.size();
    vector<char> is_bridge(m, 0), skipped(m, 0);
    for (int i = 0; i < m; i++){
        skipped[i] = 1;
        auto cut = labels(n, edges, skipped);
        is_bridge[i] = cut[edges[i].first] != cut[edges[i].second];
        skipped[i] = 0;
    }
    return is_bridge;
}

/// Bridge flags by an iterative Tarjan low-link search that skips only the edge it arrived by, O(n + m)
vector<char> tarjan_bridges(int n, const vector<pair<int, int>>& edges){
    int m = edges.size(), timer = 0;
    vector<vector<pair<int, int>>> adj(n);
    for (int i = 0; i < m; i++){
        adj[edges[i].first].push_back({edges[i].second, i});
        adj[edges[i].second].push_back({edges[i].first, i});
    }

    vector<char> is_bridge(m, 0);
    vector<int> disc(n, -1), low(n), in_edge(n, -1), it(n, 0), stack;
    for (int s = 0; s < n; s++){
        if (disc[s] != -1) continue;
        disc[s] = low[s] = timer++;
        stack.push_back(s);
        while (!stack.empty()){
            int u = stack.back();
            if (it[u] < (int)adj[u].size()){
                auto [v, id] = adj[u][it[u]++];
                if (id == in_edge[u]) continue;
                if (disc[v] == -1){
                    disc[v] = low[v] = timer++;
                    in_edge[v] = id;
                    stack.push_back(v);
                }
                else low[u] = min(low[u], disc[v]);
                continue;
            }
            stack.pop_back();
            if (stack.empty()) continue;
            int p = stack.back();
            low[p] = min(low[p], low[u]);
            if (low[u] > disc[p]) is_bridge[in_edge[u]] = 1;
        }
    }
    return is_bridge;
}

/// Bridge count, connectivity and 2-edge-connectivity of every pair against the given bridge flags
void check_all_pairs(OnlineBridges& g, int n, const vector<pair<int, int>>& edges, const vector<char>& is_bridge){
    assert(g.bridges == count(is_bridge.begin(), is_bridge.end(), 1));
    auto conn = labels(n, edges, vector<char>(edges.size(), 0));
    auto two = labels(n, edges, is_bridge);
    for (int u = 0; u < n; u++){
        int c = g.component(u);
        assert(0 <= c && c < n && two[c] == two[u]);
        for (int v = 0; v < n; v++){
            assert(g.connected(u, v) == (conn[u] == conn[v]));
            assert(g.two_edge_connected(u, v) == (two[u] == two[v]));
        }
    }
}

/// Labels must induce the same partition as the reference labels, checked through representatives in O(n)
void check_partition(OnlineBridges& g, int n, const vector<pair<int, int>>& edges, const vector<char>& is_bridge){
    assert(g.bridges == count(is_bridge.begin(), is_bridge.end(), 1));
    auto conn = labels(n, edges, vector<char>(edges.size(), 0));
    auto two = labels(n, edges, is_bridge);
    vector<int> rep_of_two(n, -1), rep_of_conn(n, -1), two_of_rep(n, -1), conn_of_cc(n, -1);
    for (int u = 0; u < n; u++){
        int c = g.component(u), cc = g.find_cc(u);
        assert(two[c] == two[u]);
        if (rep_of_two[two[u]] == -1) rep_of_two[two[u]] = c;
        assert(rep_of_two[two[u]] == c);
        if (two_of_rep[c] == -1) two_of_rep[c] = two[u];
        assert(two_of_rep[c] == two[u]);
        if (rep_of_conn[conn[u]] == -1) rep_of_conn[conn[u]] = cc;
        assert(rep_of_conn[conn[u]] == cc);
        if (conn_of_cc[cc] == -1) conn_of_cc[cc] = conn[u];
        assert(conn_of_cc[cc] == conn[u]);
    }
}

pair<int, int> random_edge(int n, int it){
    int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
    if (it % 3 == 0 && n > 1) v = (u + stress::rand_int(1, 2)) % n;
    return {u, v};
}

int main(){
    /// Small graphs, parallel edges and self loops included: every pair after every insertion against removing each edge
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 9), m = stress::rand_int(0, 3 * n);
        OnlineBridges g(n);
        vector<pair<int, int>> edges;
        for (int i = 0; i < m; i++){
            edges.push_back(random_edge(n, it));
            g.add_edge(edges.back().first, edges.back().second);
            auto is_bridge = brute_bridges(n, edges);
            assert(is_bridge == tarjan_bridges(n, edges));
            check_all_pairs(g, n, edges, is_bridge);
        }
    }

    /// Medium graphs, sparse so bridges keep appearing and vanishing: full partition after every insertion against Tarjan
    for (long long it = 0; it < stress::scaled(40); it++){
        int n = stress::rand_int(2, 400), m = stress::rand_int(n / 2, 2 * n);
        OnlineBridges g(n);
        vector<pair<int, int>> edges;
        for (int i = 0; i < m; i++){
            edges.push_back(random_edge(n, it));
            g.add_edge(edges.back().first, edges.back().second);
            check_partition(g, n, edges, tarjan_bridges(n, edges));
        }
    }

    /***
     * A path grown at alternate ends, the new vertex always passed second: it must be the re-rooted side
     * Re-rooting the growing path instead walks it end to end every time, O(n^2)
     * Then random edges merge cycles, checked against Tarjan at checkpoints
    ***/
    int n = 200000;
    OnlineBridges g(n);
    vector<pair<int, int>> edges;
    auto start = chrono::steady_clock::now();
    for (int v = 1, left = 0, right = 0; v < n; v++){
        int& end = v % 2 ? right : left;
        edges.push_back({end, v});
        g.add_edge(end, v);
        end = v;
        if (v % 4096 == 0) assert(chrono::steady_clock::now() - start < chrono::seconds(2));
    }
    assert(g.bridges == n - 1);
    for (int i = 0; i < n; i++){
        edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        g.add_edge(edges.back().first, edges.back().second);
        if (i % 50000 == 0 || i == n - 1) check_partition(g, n, edges, tarjan_bridges(n, edges));
    }

    /// Random tree edges inserted in random order, then cycles, so re-rooting walks real trees
    edges.clear();
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());
    for (int i = 1; i < n; i++) edges.push_back({order[i], order[stress::rand_int(0, i - 1)]});
    shuffle(edges.begin(), edges.end(), stress::rng());

    OnlineBridges tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);
    check_partition(tree, n, edges, tarjan_bridges(n, edges));
    for (int i = 0; i < n / 4; i++){
        edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        tree.add_edge(edges.back().first, edges.back().second);
    }
    check_partition(tree, n, edges, tarjan_bridges(n, edges));
    return 0;
}
