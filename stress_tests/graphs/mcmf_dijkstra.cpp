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

void check(int n, const vector<array<long long, 4>>& edges, long long limit){
    int s = 0, t = n - 1;
    MCMF g(n);
    vector<int> ids;
    for (auto [u, v, cap, cost] : edges) ids.push_back(g.add_edge(u, v, cap, cost));
    auto got = g.solve(s, t, limit);
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
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(2, 8), m = stress::rand_int(0, 20);
        check_slope(n, random_graph(n, m, 5, 0, 4, false));
        check_slope(n, random_graph(n, m, 5, -6, 6, true));
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(20, 60), m = stress::rand_int(n, 6 * n);
        check(n, random_graph(n, m, 1000000, 0, 1000000, false), MCMF::INF);
        check(n, random_graph(n, m, 1000000, -1000000, 1000000, true), MCMF::INF);
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
