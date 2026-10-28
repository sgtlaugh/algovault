#include "../common.h"

#define main library_main
#include "../../code_library/graphs/mcmf_spfa.cpp"
#undef main

struct Arc{ int u, v, id; long long cap, cost; };

long long min_cut(int nodes, int src, int sink, const vector<Arc>& arcs){
    long long best = LLONG_MAX;
    for (int mask = 0; mask < (1 << nodes); mask++){
        if (!(mask >> src & 1) || (mask >> sink & 1)) continue;
        long long cut = 0;
        for (auto& a : arcs) if ((mask >> a.u & 1) && !(mask >> a.v & 1)) cut += a.cap;
        best = min(best, cut);
    }

    return best;
}

/// A flow of its value is cheapest exactly when its residual graph has no negative cycle
bool has_negative_cycle(int nodes, const vector<array<long long, 3>>& residual){
    vector<long long> dist(nodes, 0);
    for (int round = 0; round <= nodes; round++){
        bool relaxed = false;
        for (auto& [u, v, w] : residual){
            if (dist[u] + w < dist[v]) dist[v] = dist[u] + w, relaxed = true;
        }
        if (!relaxed) return false;
    }

    return true;
}

/// Flow on every arc within capacity, conserved, summing to value at cost, and cheapest for its value
void check_flow(const MCMFSpfa& g, int src, int sink, const vector<Arc>& arcs, long long value, long long cost){
    long long total = 0;
    vector<long long> net(g.n, 0);
    vector<array<long long, 3>> residual;
    for (auto& a : arcs){
        long long f = g.flow(a.id);
        assert(0 <= f && f <= a.cap);
        net[a.u] -= f, net[a.v] += f, total += f * a.cost;
        if (f < a.cap) residual.push_back({a.u, a.v, a.cost});
        if (f > 0) residual.push_back({a.v, a.u, -a.cost});
    }
    for (int x = 0; x < g.n; x++) assert(net[x] == (x == sink ? value : x == src ? -value : 0));
    assert(cost == total);
    assert(!has_negative_cycle(g.n, residual));
}

int main(){
    for (long long it = 0; it < stress::scaled(15000); it++){
        int n = stress::rand_int(2, it % 10 ? 8 : 12);
        int src = stress::rand_int(0, n - 1), sink = (src + stress::rand_int(1, n - 1)) % n;

        /// cost = base + p[u] - p[v] with base >= 0 allows negative arcs but no negative cycle
        vector<long long> potential(n);
        for (auto& p : potential) p = it % 2 ? stress::rand_int(0, 20) : 0;

        MCMFSpfa g(n);
        vector<Arc> arcs;
        for (int m = stress::rand_int(0, 2 * n + 4); m; m--){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long cap = stress::rand_int(0, it % 5 ? 8 : 1000000000), base = stress::rand_int(0, 15);
            if (stress::rand_int(0, 3)){
                long long cost = base + potential[u] - potential[v];
                arcs.push_back({u, v, g.add_edge(u, v, cap, cost), cap, cost});
            }
            else{
                long long cost = base + abs(potential[u] - potential[v]);  /// reduced cost stays >= 0 both ways
                arcs.push_back({u, v, g.add_edge(u, v, cap, cost), cap, cost});
                arcs.push_back({v, u, g.add_edge(v, u, cap, cost), cap, cost});
            }
        }
        MCMFSpfa limited = g;

        auto [flow, cost] = g.solve(src, sink);
        assert(flow == min_cut(n, src, sink, arcs));
        check_flow(g, src, sink, arcs, flow, cost);

        /// A limited solve stops at the cheapest flow of that value, a second solve continues to the full answer
        long long k = stress::rand_int(0, flow);
        auto [first_flow, first_cost] = limited.solve(src, sink, k);
        assert(first_flow == k);
        check_flow(limited, src, sink, arcs, k, first_cost);
        auto [rest_flow, rest_cost] = limited.solve(src, sink);
        assert(first_flow + rest_flow == flow && first_cost + rest_cost == cost);
    }

    return 0;
}
