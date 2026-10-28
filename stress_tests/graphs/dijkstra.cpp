#include "../common.h"

#define main library_main
#include "../../code_library/graphs/dijkstra.cpp"
#undef main

struct TestEdge{
    int u, v;
    long long w;
    bool directed;
};

/// lightest[{u, v}] is the weight of the lightest u -> v arc, an undirected edge adds both arcs
map<pair<int, int>, long long> lightest_arcs(const vector<TestEdge>& edges){
    map<pair<int, int>, long long> lightest;
    auto add = [&](int u, int v, long long w){
        auto it = lightest.find({u, v});
        if (it == lightest.end() || w < it->second) lightest[{u, v}] = w;
    };

    for (auto e : edges){
        add(e.u, e.v, e.w);
        if (!e.directed) add(e.v, e.u, e.w);
    }

    return lightest;
}

/// The oracles keep their own sentinel so a library DIJKSTRA_INF that is too small cannot cap the expected distances
const long long UNREACHED = numeric_limits<long long>::max();

vector<long long> to_library_inf(vector<long long> d){
    for (auto& x : d) if (x == UNREACHED) x = DIJKSTRA_INF;
    return d;
}

vector<vector<long long>> floyd(int n, const map<pair<int, int>, long long>& lightest){
    vector<vector<long long>> d(n, vector<long long>(n, UNREACHED));
    for (int v = 0; v < n; v++) d[v][v] = 0;
    for (auto [arc, w] : lightest) d[arc.first][arc.second] = min(d[arc.first][arc.second], w);

    for (int k = 0; k < n; k++){
        for (int i = 0; i < n; i++){
            if (d[i][k] == UNREACHED) continue;
            for (int j = 0; j < n; j++){
                if (d[k][j] != UNREACHED) d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
            }
        }
    }
    for (auto& row : d) row = to_library_inf(row);
    return d;
}

vector<long long> bellman_ford(int n, const map<pair<int, int>, long long>& lightest, int src){
    vector<long long> d(n, UNREACHED);
    d[src] = 0;
    for (bool changed = true; changed; ){
        changed = false;
        for (auto [arc, w] : lightest){
            auto [u, v] = arc;
            if (d[u] != UNREACHED && d[u] + w < d[v]) d[v] = d[u] + w, changed = true;
        }
    }

    return to_library_inf(d);
}

/// One vertex at a time, since after run(src, target) only v = target is promised
template<typename Graph>
void check_vertex(const Graph& g, const map<pair<int, int>, long long>& lightest, const vector<int>& sources, int v,
                  long long expected){
    assert(g.dist[v] == expected);
    vector<int> path = g.path(v);
    if (expected == DIJKSTRA_INF){
        assert(path.empty());
        return;
    }

    assert(!path.empty() && count(sources.begin(), sources.end(), path.front()) && path.back() == v);
    long long total = 0;
    for (size_t i = 0; i + 1 < path.size(); i++){
        auto it = lightest.find({path[i], path[i + 1]});
        assert(it != lightest.end());
        total += it->second;
    }
    assert(total == expected);
}

template<typename Graph>
void check_run(const Graph& g, const map<pair<int, int>, long long>& lightest, const vector<int>& sources,
               const vector<long long>& expected){
    assert(g.dist == expected);
    for (int v = 0; v < g.n; v++) check_vertex(g, lightest, sources, v, expected[v]);
}

/// A few target runs first, then a full run that must not see their leftovers
template<typename Graph>
void check_targets_then_full(Graph& g, const map<pair<int, int>, long long>& lightest, int src, const vector<long long>& expected){
    for (int i = 0; i < 5; i++){
        int target = stress::rand_int(0, g.n - 1);
        g.run(src, target);
        check_vertex(g, lightest, {src}, target, expected[target]);
    }
    g.run(src);
    check_run(g, lightest, {src}, expected);
}

template<typename Graph>
Graph build(int n, const vector<TestEdge>& edges){
    Graph g(n);
    for (auto e : edges) g.add_edge(e.u, e.v, e.w, e.directed);
    return g;
}

/// Every source on one graph, target runs between full runs, so a full run after an early exit must reset everything
template<typename Graph>
void check_all_sources(int n, const vector<TestEdge>& edges){
    Graph g = build<Graph>(n, edges);
    auto lightest = lightest_arcs(edges);
    auto expected = floyd(n, lightest);
    for (int src = 0; src < n; src++){
        for (int target = 0; target < n; target++){
            g.run(src, target);
            check_vertex(g, lightest, {src}, target, expected[src][target]);
        }
        g.run(src);
        check_run(g, lightest, {src}, expected[src]);
    }

    /// Repeated sources included, the distance to v is the minimum over the sources
    for (int i = 0; i < 3; i++){
        vector<int> sources(stress::rand_int(1, n));
        for (int& src : sources) src = stress::rand_int(0, n - 1);
        vector<long long> nearest(n, DIJKSTRA_INF);
        for (int src : sources){
            for (int v = 0; v < n; v++) nearest[v] = min(nearest[v], expected[src][v]);
        }

        int target = stress::rand_int(0, n - 1);
        g.run(sources, target);
        check_vertex(g, lightest, sources, target, nearest[target]);
        g.run(sources);
        check_run(g, lightest, sources, nearest);
    }
}

/// Every arc as a callback, checked against the explicit oracle, with neighbors called at most once per state
vector<long long> run_implicit(int n, const vector<TestEdge>& edges, const vector<int>& sources, int target = -1){
    vector<vector<pair<int, long long>>> adj(n);
    for (auto e : edges){
        adj[e.u].push_back({e.v, e.w});
        if (!e.directed) adj[e.v].push_back({e.u, e.w});
    }

    vector<long long> seen(n, -1);
    auto dist = dijkstra_implicit(n, sources, [&](int u, long long du, DijkstraRelax& relax){
        assert(seen[u] == -1);
        seen[u] = du;
        for (auto [v, w] : adj[u]) relax(v, du + w);
    }, target);
    for (int v = 0; v < n; v++) assert(seen[v] == -1 || seen[v] == dist[v]);
    return dist;
}

void check_implicit_all_sources(int n, const vector<TestEdge>& edges){
    auto expected = floyd(n, lightest_arcs(edges));
    for (int src = 0; src < n; src++){
        assert(run_implicit(n, edges, {src}) == expected[src]);
        int target = stress::rand_int(0, n - 1);
        assert(run_implicit(n, edges, {src}, target)[target] == expected[src][target]);
    }

    vector<int> sources(stress::rand_int(1, n));
    for (int& src : sources) src = stress::rand_int(0, n - 1);
    vector<long long> nearest(n, DIJKSTRA_INF);
    for (int src : sources){
        for (int v = 0; v < n; v++) nearest[v] = min(nearest[v], expected[src][v]);
    }
    assert(run_implicit(n, edges, sources) == nearest);
}

/// A run after more add_edge calls must see the new edges, which Dijkstra's cached adjacency could miss
template<typename Graph>
void incremental_case(int n, const vector<TestEdge>& edges){
    Graph g(n);
    vector<TestEdge> added;
    for (auto e : edges){
        g.add_edge(e.u, e.v, e.w, e.directed);
        added.push_back(e);
        if (stress::rand_int(0, 2) > 0) continue;

        auto lightest = lightest_arcs(added);
        int src = stress::rand_int(0, n - 1);
        g.run(src);
        check_run(g, lightest, {src}, floyd(n, lightest)[src]);
    }
}

/// On the chain 0 -> 1 -> ... -> n - 1 a search that stops when it settles k has never relaxed the arc k -> k + 1
template<typename Graph>
void early_exit_case(){
    int n = stress::rand_int(2, 50);
    Graph g(n);
    for (int v = 0; v + 1 < n; v++) g.add_edge(v, v + 1, stress::rand_int(0, 1));

    int target = stress::rand_int(0, n - 2);
    g.run(0, target);
    assert(g.dist[target + 1] == DIJKSTRA_INF);

    g.run(0);
    assert(g.dist[n - 1] != DIJKSTRA_INF);
}

/// 0 fans out to k vertices before its 0 edge to target, so the search stops with the fan still queued, unexpanded
template<typename Graph>
void wide_early_exit_case(){
    int k = stress::rand_int(1, 30), target = k + 1;
    Graph g(2 * k + 2);
    for (int i = 1; i <= k; i++){
        g.add_edge(0, i, 1);
        g.add_edge(i, target + i, 0);
    }
    g.add_edge(0, target, 0);

    g.run(0, target);
    assert(g.dist[target] == 0);
    for (int i = 1; i <= k; i++) assert(g.dist[target + i] == DIJKSTRA_INF);
}

/// An 11 edge chain near 2e17 per step, shortcuts at most 2e17 per skipped step: distances reach, never pass, 2.2e18
long long heavy_chain_case(){
    int n = 12;
    const long long step = 200000000000000000LL;
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());

    vector<TestEdge> edges;
    for (int i = 0; i + 1 < n; i++){
        long long w = stress::rand_int(step - step / 20, step);
        edges.push_back({order[i], order[i + 1], w, stress::rand_int(0, 1) == 0});
    }
    for (int s = stress::rand_int(0, n); s > 0; s--){
        int i = stress::rand_int(0, n - 2), j = stress::rand_int(i + 1, n - 1);
        long long w = stress::rand_int((j - i) * (step - step / 20), (j - i) * step);
        edges.push_back({order[i], order[j], w, true});
    }

    check_all_sources<Dijkstra>(n, edges);
    check_all_sources<DenseDijkstra>(n, edges);
    check_implicit_all_sources(n, edges);

    long long longest = 0;
    for (auto& row : floyd(n, lightest_arcs(edges))){
        for (long long d : row) if (d != DIJKSTRA_INF) longest = max(longest, d);
    }

    return longest;
}

/// The 1 chain a_0 -> ... -> a_k feeds a_i into the 0 chain at z_{2i}, so a_i's wave reaches z_p before a_{i-1}'s
void zero_one_adversarial_case(){
    int k = 20000, n = 3 * k + 1;
    auto a = [&](int i){ return i; };
    auto z = [&](int p){ return k + p; };

    vector<TestEdge> edges;
    for (int i = 0; i < k; i++) edges.push_back({a(i), a(i + 1), 1, true});
    for (int i = 1; i <= k; i++) edges.push_back({a(i), z(2 * i), 0, true});
    for (int p = 1; p < 2 * k; p++) edges.push_back({z(p), z(p + 1), 0, true});

    auto heap = build<Dijkstra>(n, edges);
    auto zero_one = build<ZeroOneBFS>(n, edges);
    heap.run(0);

    /// a FIFO queue relaxes z_p about p / 2 times, Theta(k^2) total, a deque O(k)
    auto start = chrono::steady_clock::now();
    zero_one.run(0);
    assert(chrono::steady_clock::now() - start < chrono::seconds(2));
    assert(zero_one.dist == heap.dist);
}

/// The implicit twin of early_exit_case: settling target must stop before expanding it
void implicit_early_exit_case(){
    for (int n = 2; n <= 50; n++){
        vector<TestEdge> chain;
        for (int v = 0; v + 1 < n; v++) chain.push_back({v, v + 1, stress::rand_int(0, 1), true});

        int target = stress::rand_int(0, n - 2);
        assert(run_implicit(n, chain, {0}, target)[target + 1] == DIJKSTRA_INF);
        assert(run_implicit(n, chain, {0})[n - 1] != DIJKSTRA_INF);
    }
}

/// Arc u -> v turns du into du + add + du * pct / 100, growing with du, so the best route depends on the load so far
struct TollArc{
    int u, v;
    long long add, pct;

    long long apply(long long du) const{
        return du + add + du * pct / 100;
    }
};

/// Since every arc only adds, a best walk is a simple path, so trying every simple path from every source is exact
void simple_paths(int u, long long du, int on_path, const vector<TollArc>& arcs, vector<long long>& best){
    best[u] = min(best[u], du);
    for (auto& a : arcs){
        if (a.u == u && !(on_path >> a.v & 1)) simple_paths(a.v, a.apply(du), on_path | 1 << a.v, arcs, best);
    }
}

void toll_case(){
    int n = stress::rand_int(1, 7), m = stress::rand_int(0, 4 * n);
    vector<TollArc> arcs(m);
    for (auto& a : arcs){
        a.u = stress::rand_int(0, n - 1), a.v = stress::rand_int(0, n - 1);
        a.add = stress::rand_int(0, 3) == 0 ? 0 : stress::rand_int(0, 30);
        a.pct = stress::rand_int(0, 2) == 0 ? 0 : stress::rand_int(0, 150);
    }
    vector<int> sources(stress::rand_int(1, 2));
    for (int& src : sources) src = stress::rand_int(0, n - 1);

    vector<long long> best(n, DIJKSTRA_INF);
    for (int src : sources) simple_paths(src, 0, 1 << src, arcs, best);

    auto toll = [&](int u, long long du, DijkstraRelax& relax){
        for (auto& a : arcs) if (a.u == u) relax(a.v, a.apply(du));
    };
    assert(dijkstra_implicit(n, sources, toll) == best);
    int target = stress::rand_int(0, n - 1);
    assert(dijkstra_implicit(n, sources, toll, target)[target] == best[target]);
}

vector<TestEdge> random_edges(int n, int m, long long max_w){
    vector<TestEdge> edges;
    for (int i = 0; i < m; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        long long w = stress::rand_int(0, 3) == 0 ? 0 : stress::rand_int(0, max_w);
        edges.push_back({u, v, w, stress::rand_int(0, 1) == 0});
    }

    return edges;
}

int main(){
    /// 2e17 keeps every distance on at most 11 edges below 2.2e18 < DIJKSTRA_INF, the documented limit
    const long long weight_limits[] = {1, 5, 1000000, 200000000000000000LL};
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 12), m = stress::rand_int(0, 3 * n);
        long long max_w = weight_limits[it % 4];
        auto edges = random_edges(n, m, max_w);

        check_all_sources<Dijkstra>(n, edges);
        check_all_sources<DenseDijkstra>(n, edges);
        incremental_case<Dijkstra>(n, edges);
        incremental_case<DenseDijkstra>(n, edges);
        if (max_w > 1) continue;
        check_all_sources<ZeroOneBFS>(n, edges);
        incremental_case<ZeroOneBFS>(n, edges);
    }

    for (long long it = 0; it < stress::scaled(12); it++){
        int n = stress::rand_int(500, 2000), m = stress::rand_int(n, 6 * n);
        long long max_w = it % 2 ? 1 : 1000000000;
        auto edges = random_edges(n, m, max_w);
        auto lightest = lightest_arcs(edges);
        int src = stress::rand_int(0, n - 1);
        auto expected = bellman_ford(n, lightest, src);

        auto sparse = build<Dijkstra>(n, edges);
        check_targets_then_full(sparse, lightest, src, expected);

        auto dense = build<DenseDijkstra>(n, edges);
        check_targets_then_full(dense, lightest, src, expected);

        int target = (src + n / 2) % n;
        assert(run_implicit(n, edges, {src}, target)[target] == expected[target]);
        assert(run_implicit(n, edges, {src}) == sparse.dist);

        if (max_w > 1) continue;
        auto zero_one = build<ZeroOneBFS>(n, edges);
        check_targets_then_full(zero_one, lightest, src, expected);
    }

    /// A 300 x 300 grid of 0/1 weights, 0-1 BFS against Dijkstra on 90000 vertices
    int side = 300;
    vector<TestEdge> grid;
    for (int r = 0; r < side; r++){
        for (int c = 0; c < side; c++){
            int v = r * side + c;
            if (c + 1 < side) grid.push_back({v, v + 1, stress::rand_int(0, 1), false});
            if (r + 1 < side) grid.push_back({v, v + side, stress::rand_int(0, 1), false});
        }
    }

    auto heap_grid = build<Dijkstra>(side * side, grid);
    auto deque_grid = build<ZeroOneBFS>(side * side, grid);
    heap_grid.run(0), deque_grid.run(0);
    assert(heap_grid.dist == deque_grid.dist);
    assert(run_implicit(side * side, grid, {0}) == heap_grid.dist);

    long long longest = 0;
    for (long long it = 0; it < stress::scaled(300); it++) longest = max(longest, heavy_chain_case());
    assert(longest > 2000000000000000000LL);

    zero_one_adversarial_case();
    for (long long it = 0; it < stress::scaled(200); it++){
        early_exit_case<Dijkstra>();
        early_exit_case<DenseDijkstra>();
        early_exit_case<ZeroOneBFS>();
        wide_early_exit_case<Dijkstra>();
        wide_early_exit_case<DenseDijkstra>();
        wide_early_exit_case<ZeroOneBFS>();
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 12);
        check_implicit_all_sources(n, random_edges(n, stress::rand_int(0, 3 * n), weight_limits[it % 4]));
    }
    for (long long it = 0; it < stress::scaled(5000); it++) toll_case();
    implicit_early_exit_case();
    /// A dv below du is a negative arc, here at du = 5 so an assert against a stale du of 0 would miss it
    stress::expect_assert("graphs/dijkstra.cpp", []{
        dijkstra_implicit(3, {0}, [](int u, long long du, DijkstraRelax& relax){
            if (u == 0) relax(1, 5);
            if (u == 1) relax(2, du - 1);
        });
    });
    return 0;
}
