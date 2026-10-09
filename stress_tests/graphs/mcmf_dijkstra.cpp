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

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(2, 8), m = stress::rand_int(0, 20);
        check(n, random_graph(n, m, 5, 0, 10, false), MCMF::INF);
        check(n, random_graph(n, m, 5, -10, 10, true), MCMF::INF);
        check(n, random_graph(n, m, 5, 0, 10, false), stress::rand_int(0, 6));
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
