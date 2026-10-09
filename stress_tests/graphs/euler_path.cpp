#include "../common.h"

#define main library_main
#include "../../code_library/graphs/euler_path.cpp"
#undef main

struct Edge{
    int u, v;
};

/// The start the header promises for start = -1: first out - in = 1 or odd degree vertex, else first vertex with an edge
template<bool directed>
int default_start(int n, const vector<Edge>& edges){
    vector<int> out(n, 0), in(n, 0);
    for (auto [u, v] : edges) out[u]++, in[v]++;
    for (int v = 0; v < n; v++){
        if (directed ? out[v] - in[v] == 1 : (out[v] + in[v]) % 2 == 1) return v;
    }
    for (int v = 0; v < n; v++){
        if (out[v] + in[v] > 0) return v;
    }
    return 0;
}

/// Every edge exactly once, consecutive edges share endpoints, the walk begins at start (or the promised default)
template<bool directed>
void check_walk(int n, const vector<Edge>& edges, int start, const vector<int>& path, const vector<int>& walk){
    int m = edges.size();
    assert((int)path.size() == m + 1 && (int)walk.size() == m);
    assert(path[0] == (start == -1 ? default_start<directed>(n, edges) : start));
    for (int v : path) assert(0 <= v && v < n);

    vector<char> seen(m, 0);
    for (int i = 0; i < m; i++){
        int id = walk[i];
        assert(0 <= id && id < m && !seen[id]);
        seen[id] = 1;
        auto [u, v] = edges[id];
        bool forward = u == path[i] && v == path[i + 1];
        bool backward = v == path[i] && u == path[i + 1];
        assert(forward || (!directed && backward));
    }
}

/// Exhaustive search over edge orderings: can a trail from start use every edge?
template<bool directed>
bool brute_from(int u, const vector<Edge>& edges, vector<char>& used, int left){
    if (left == 0) return true;
    for (int id = 0; id < (int)edges.size(); id++){
        if (used[id]) continue;
        auto [a, b] = edges[id];
        int next = a == u ? b : (!directed && b == u ? a : -1);
        if (next == -1) continue;
        used[id] = 1;
        bool ok = brute_from<directed>(next, edges, used, left - 1);
        used[id] = 0;
        if (ok) return true;
    }
    return false;
}

/// The textbook condition: degrees balanced up to one source/sink pair, every vertex with an edge in one weak component
template<bool directed>
bool condition_holds(int n, const vector<Edge>& edges, int start){
    vector<int> parent(n), out(n, 0), in(n, 0);
    iota(parent.begin(), parent.end(), 0);
    function<int(int)> find = [&](int x){ return parent[x] == x ? x : parent[x] = find(parent[x]); };
    for (auto [u, v] : edges) out[u]++, in[v]++, parent[find(u)] = find(v);

    int root = -1, plus = 0, minus = 0, odd = 0;
    for (int v = 0; v < n; v++){
        if (out[v] + in[v] == 0) continue;
        if (root == -1) root = find(v);
        if (find(v) != root) return false;
        if (directed){
            int b = out[v] - in[v];
            if (b > 1 || b < -1) return false;
            plus += b == 1, minus += b == -1;
        }
        else odd += (out[v] + in[v]) % 2;
    }
    if (edges.empty()) return start < n;

    if (directed){
        if (plus != minus || plus > 1) return false;
        if (start == -1) return true;
        return plus ? out[start] - in[start] == 1 : out[start] > 0;
    }
    if (odd != 0 && odd != 2) return false;
    if (start == -1) return true;
    return odd ? (out[start] + in[start]) % 2 == 1 : out[start] + in[start] > 0;
}

vector<Edge> random_edges(int n, int m){
    vector<Edge> edges(m);
    for (auto& [u, v] : edges) u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
    return edges;
}

template<bool directed>
void check(int n, const vector<Edge>& edges, bool brute){
    EulerPath<directed> g(n);
    for (int i = 0; i < (int)edges.size(); i++) assert(g.add_edge(edges[i].u, edges[i].v) == i);

    bool any = false;
    for (int start = n - 1; start >= -1; start--){
        auto path = g.solve(start);
        bool expected = condition_holds<directed>(n, edges, start);
        if (brute && start != -1){
            vector<char> used(edges.size(), 0);
            assert(brute_from<directed>(start, edges, used, edges.size()) == expected);
            any |= expected;
        }
        if (brute && start == -1) assert(expected == any);
        assert(!path.empty() == expected);
        if (expected) check_walk<directed>(n, edges, start, path, g.walk_edges);
        else assert(g.walk_edges.empty());
    }
}

/// A random walk is always an Euler path of the edges it traversed
vector<Edge> random_walk(int n, int m){
    vector<Edge> edges;
    int u = stress::rand_int(0, n - 1);
    for (int i = 0; i < m; i++){
        int v = stress::rand_int(0, n - 1);
        edges.push_back({u, v});
        u = v;
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

int main(){
    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, 5), m = stress::rand_int(0, 7);
        auto edges = it % 2 ? random_walk(n, m) : random_edges(n, m);
        check<true>(n, edges, true);
        check<false>(n, edges, true);
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 40), m = stress::rand_int(0, 80);
        auto edges = it % 3 ? random_walk(n, m) : random_edges(n, m);
        check<true>(n, edges, false);
        check<false>(n, edges, false);
    }

    /// Random walks with a few extra isolated vertices, plus one long cycle and one long path for stack depth
    for (int n : {1000, 300000}){
        auto edges = random_walk(n, 2 * n);
        for (auto& e : edges) e.u = e.u * 2 % (n + 7), e.v = e.v * 2 % (n + 7);
        int total = n + 7;
        EulerDirected d(total);
        EulerUndirected u(total);
        for (auto [a, b] : edges) d.add_edge(a, b), u.add_edge(a, b);
        check_walk<true>(total, edges, -1, d.solve(), d.walk_edges);
        check_walk<false>(total, edges, -1, u.solve(), u.walk_edges);

        vector<Edge> cycle, path;
        for (int i = 0; i < n; i++) cycle.push_back({i, (i + 1) % n});
        for (int i = 0; i + 1 < n; i++) path.push_back({i, i + 1});
        EulerDirected dc(n), dp(n);
        EulerUndirected uc(n), up(n);
        for (auto [a, b] : cycle) dc.add_edge(a, b), uc.add_edge(a, b);
        for (auto [a, b] : path) dp.add_edge(a, b), up.add_edge(a, b);
        check_walk<true>(n, cycle, 0, dc.solve(), dc.walk_edges);
        check_walk<false>(n, cycle, n - 1, uc.solve(n - 1), uc.walk_edges);
        check_walk<true>(n, path, 0, dp.solve(), dp.walk_edges);
        check_walk<false>(n, path, n - 1, up.solve(n - 1), up.walk_edges);
        assert(dp.solve(n - 1).empty());
    }
    return 0;
}
