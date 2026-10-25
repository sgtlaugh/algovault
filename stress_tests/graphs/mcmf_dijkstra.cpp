#include "../common.h"

#define main library_main
#include "../../code_library/graphs/mcmf_dijkstra.cpp"
#undef main

/// Successive shortest paths with a full Bellman Ford each round, no potentials
pair<long long, long long> brute(int n, const vector<array<long long, 4>>& input, int s, int t, long long limit){
    struct E{ int to; long long cap, cost; };
    vector<E> e;
    vector<vector<int>> adj(n);
    for (auto [u, v, cap, cost] : input){
        adj[u].push_back(e.size()), e.push_back({(int)v, cap, cost});
        adj[v].push_back(e.size()), e.push_back({(int)u, 0, -cost});
    }

    long long flow = 0, cost = 0;
    while (flow < limit){
        vector<long long> dist(n, LLONG_MAX);
        vector<int> par(n, -1);
        dist[s] = 0;
        for (int round = 0; round < n; round++){
            for (int u = 0; u < n; u++){
                if (dist[u] == LLONG_MAX) continue;
                for (int id : adj[u]){
                    if (e[id].cap > 0 && dist[u] + e[id].cost < dist[e[id].to]) dist[e[id].to] = dist[u] + e[id].cost, par[e[id].to] = id;
                }
            }
        }

        if (dist[t] == LLONG_MAX) break;
        long long push = limit - flow;
        for (int v = t; v != s; v = e[par[v] ^ 1].to) push = min(push, e[par[v]].cap);
        for (int v = t; v != s; v = e[par[v] ^ 1].to) e[par[v]].cap -= push, e[par[v] ^ 1].cap += push, cost += push * e[par[v]].cost;
        flow += push;
    }
    return {flow, cost};
}

vector<array<long long, 4>> random_graph(int n, int m, long long max_cap, long long min_cost, long long max_cost, bool dag){
    /// DAG labels are shuffled (source and sink stay put) so Bellman Ford cannot converge in one index-order pass
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    if (n > 2) shuffle(label.begin() + 1, label.end() - 1, stress::rng());

    vector<array<long long, 4>> edges;
    for (int i = 0; i < m; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (u == v) continue;
        if (dag && u > v) swap(u, v);
        edges.push_back({label[u], label[v], stress::rand_int(0, max_cap), stress::rand_int(min_cost, max_cost)});
    }
    return edges;
}

/// A path numbered against index order, two parallel edges per step, forward shortcuts so negative costs make no cycle
vector<array<long long, 4>> reversed_chain(int n, int extra, long long min_cost, long long max_cost){
    vector<int> path = {0};
    for (int v = n - 2; v >= 1; v--) path.push_back(v);
    path.push_back(n - 1);

    vector<array<long long, 4>> edges;
    for (int i = 0; i + 1 < n; i++){
        long long d = stress::rand_int(min_cost, max_cost);
        edges.push_back({path[i], path[i + 1], 1, d});
        edges.push_back({path[i], path[i + 1], stress::rand_int(0, 2), d + stress::rand_int(0, max_cost - min_cost)});
    }
    for (int i = 0; i < extra; i++){
        int a = stress::rand_int(0, n - 2), b = stress::rand_int(a + 1, n - 1);
        edges.push_back({path[a], path[b], stress::rand_int(0, 3), stress::rand_int(min_cost, max_cost)});
    }
    return edges;
}

/// A negative reduced cost lets Dijkstra settle the sink too early, so each residual edge reachable from s needs >= 0
void check_reduced_costs(const MCMF& g, int s){
    vector<char> seen(g.n, 0);
    vector<int> stack = {s};
    seen[s] = 1;
    while (!stack.empty()){
        int u = stack.back();
        stack.pop_back();
        assert(g.potential[u] < MCMF::INF);
        for (int id : g.adj[u]){
            const auto& e = g.edges[id];
            if (e.cap <= 0) continue;
            assert(g.potential[u] + e.cost - g.potential[e.to] >= 0);
            if (!seen[e.to]) seen[e.to] = 1, stack.push_back(e.to);
        }
    }
}

void check(int n, const vector<array<long long, 4>>& edges, long long limit){
    int s = 0, t = n - 1;
    MCMF g(n);
    vector<int> ids;
    for (auto [u, v, cap, cost] : edges) ids.push_back(g.add_edge(u, v, cap, cost));
    g.init_potential(s);
    check_reduced_costs(g, s);
    auto got = g.solve(s, t, limit);
    MCMF residual = g;
    residual.init_potential(s);
    check_reduced_costs(residual, s);
    assert(got == brute(n, edges, s, t, limit));

    vector<long long> balance(n, 0);
    long long cost = 0;
    for (int i = 0; i < (int)edges.size(); i++){
        long long f = g.flow(ids[i]);
        assert(0 <= f && f <= edges[i][2]);
        balance[edges[i][0]] -= f, balance[edges[i][1]] += f, cost += f * edges[i][3];
    }
    for (int v = 0; v < n; v++) assert(balance[v] == (v == t ? got.first : v == s ? -got.first : 0));
    assert(cost == got.second);
}

/// Each round stops at the sink, so its potential update must still leave reduced costs >= 0 and the found path tight
void check_rounds(int n, const vector<array<long long, 4>>& edges){
    int s = 0, t = n - 1;
    MCMF g(n);
    for (auto [u, v, cap, cost] : edges) g.add_edge(u, v, cap, cost);
    g.init_potential(s);

    long long flow = 0, cost = 0;
    while (g.dijkstra(s, t)){
        check_reduced_costs(g, s);
        long long push = MCMF::INF;
        for (int v = t; v != s; v = g.edges[g.parent_edge[v] ^ 1].to){
            int id = g.parent_edge[v], u = g.edges[id ^ 1].to;
            assert(g.potential[u] + g.edges[id].cost - g.potential[v] == 0);
            push = min(push, g.edges[id].cap);
        }
        for (int v = t; v != s; v = g.edges[g.parent_edge[v] ^ 1].to){
            int id = g.parent_edge[v];
            g.edges[id].cap -= push, g.edges[id ^ 1].cap += push, cost += push * g.edges[id].cost;
        }
        flow += push;
        check_reduced_costs(g, s);
    }

    assert((make_pair(flow, cost) == brute(n, edges, s, t, MCMF::INF)));
}

/// Breakpoints against a fresh brute force solve at every flow limit from 0 to past the max flow
void check_slope(int n, const vector<array<long long, 4>>& edges){
    int s = 0, t = n - 1;
    MCMF g(n);
    for (auto [u, v, cap, cost] : edges) g.add_edge(u, v, cap, cost);
    auto points = g.slope(s, t);

    assert((points[0] == make_pair(0LL, 0LL)));
    for (int i = 1; i < (int)points.size(); i++){
        long long dx = points[i].first - points[i - 1].first, dy = points[i].second - points[i - 1].second;
        assert(dx > 0 && dy % dx == 0);
        if (i > 1) assert(dy / dx > (points[i - 1].second - points[i - 2].second) / (points[i - 1].first - points[i - 2].first));
    }

    long long max_flow = points.back().first;
    assert(brute(n, edges, s, t, MCMF::INF).first == max_flow);
    for (long long limit = 0, i = 0; limit <= max_flow + 1; limit++){
        while (i + 1 < (int)points.size() && points[i + 1].first < limit) i++;
        auto expected = brute(n, edges, s, t, limit);
        long long f = min(limit, max_flow), cost = points[i].second;
        if (f > points[i].first){
            auto [x0, y0] = points[i];
            auto [x1, y1] = points[i + 1];
            cost += (f - x0) * ((y1 - y0) / (x1 - x0));
        }
        assert((make_pair(f, cost) == expected));

        MCMF h(n);
        for (auto [u, v, cap, c] : edges) h.add_edge(u, v, cap, c);
        auto truncated = h.slope(s, t, limit);
        assert(truncated.back() == expected);
        for (int j = 0; j + 1 < (int)truncated.size(); j++) assert(truncated[j] == points[j]);

        auto rest = h.solve(s, t);
        assert((make_pair(expected.first + rest.first, expected.second + rest.second) == points.back()));
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(2, 8), m = stress::rand_int(0, 20);
        check(n, random_graph(n, m, 5, 0, 10, false), MCMF::INF);
        check(n, random_graph(n, m, 5, -10, 10, true), MCMF::INF);
        check(n, random_graph(n, m, 5, 0, 10, false), stress::rand_int(0, 6));
        check(n, random_graph(n, m, 5, -1, 3, true), stress::rand_int(0, 6));
        check_rounds(n, random_graph(n, m, 5, 0, 10, false));
        check_rounds(n, random_graph(n, m, 5, -10, 10, true));
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(2, 8), m = stress::rand_int(0, 20);
        check_slope(n, random_graph(n, m, 5, 0, 4, false));
        check_slope(n, random_graph(n, m, 5, -6, 6, true));
        check_slope(n, random_graph(n, m, 5, -1, 3, true));
        check_slope(n, reversed_chain(n, m, -6, 6));
    }

    for (long long it = 0; it < stress::scaled(500); it++){
        int n = stress::rand_int(2, 40), extra = stress::rand_int(0, 2 * n);
        check(n, reversed_chain(n, extra, 0, 10), stress::rand_int(1, 3));
        check(n, reversed_chain(n, extra, -10, 10), MCMF::INF);
        check_rounds(n, reversed_chain(n, extra, -10, 10));
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(20, 60), m = stress::rand_int(n, 6 * n);
        check(n, random_graph(n, m, 1000000, 0, 1000000, false), MCMF::INF);
        check(n, random_graph(n, m, 1000000, -1000000, 1000000, true), MCMF::INF);
        check_rounds(n, random_graph(n, m, 1000000, 0, 1000000, false));
        check_rounds(n, random_graph(n, m, 1000000, -1000000, 1000000, true));
    }

    /// Large negative-cost DAG: valid starting potentials keep this fast, no brute force at this size
    for (long long it = 0; it < stress::scaled(2); it++){
        int n = 3000;
        auto edges = random_graph(n, 30000, 1000, -1000, 1000, true);
        MCMF g(n);
        vector<int> ids;
        for (auto [u, v, cap, cost] : edges) ids.push_back(g.add_edge(u, v, cap, cost));
        auto [flow, cost] = g.solve(0, n - 1);

        long long check_cost = 0;
        vector<long long> balance(n, 0);
        for (int i = 0; i < (int)edges.size(); i++){
            long long f = g.flow(ids[i]);
            assert(0 <= f && f <= edges[i][2]);
            balance[edges[i][0]] -= f, balance[edges[i][1]] += f, check_cost += f * edges[i][3];
        }
        for (int v = 1; v + 1 < n; v++) assert(balance[v] == 0);
        assert(balance[n - 1] == flow && check_cost == cost);
    }
    return 0;
}
