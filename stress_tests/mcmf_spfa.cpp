#include "common.h"

#define main library_main
#include "../code_library/mcmf_spfa.cpp"
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

int main(){
    auto* g = new FlowGraph<long long>(1);  /// fixed size arrays, too big for the stack, rebuilt in place per case

    for (long long it = 0; it < stress::scaled(15000); it++){
        int n = stress::rand_int(1, it % 10 ? 7 : 11), nodes = n + 1;  /// the library sizes for nodes 0 .. n, so both 0 and 1 based work
        int src = stress::rand_int(0, n), sink = (src + stress::rand_int(1, n)) % nodes;

        /// cost = base + p[u] - p[v] with base >= 0 allows negative arcs but no negative cycle
        vector<long long> potential(nodes);
        for (auto& p : potential) p = it % 2 ? stress::rand_int(0, 20) : 0;

        g->~FlowGraph();
        g = new (g) FlowGraph<long long>(n);  /// the returned pointer is valid for the new object even with a const member
        vector<Arc> arcs;
        int next_id = 2;
        for (int m = stress::rand_int(0, 2 * n + 6); m; m--){
            int u = stress::rand_int(0, n), v = stress::rand_int(0, n);
            long long cap = stress::rand_int(0, it % 5 ? 8 : 1000000000), base = stress::rand_int(0, 15);
            if (stress::rand_int(0, 3)){
                long long cost = base + potential[u] - potential[v];
                g->add_edge(u, v, cap, cost);
                arcs.push_back({u, v, next_id, cap, cost}), next_id += 2;
            }
            else{
                long long cost = base + abs(potential[u] - potential[v]);  /// reduced cost stays >= 0 both ways
                g->add_edge(u, v, cap, cost, false);
                arcs.push_back({u, v, next_id, cap, cost}), arcs.push_back({v, u, next_id + 2, cap, cost}), next_id += 4;
            }
        }

        auto [cost, flow] = g->mincost_maxflow(src, sink);
        assert(flow == min_cut(nodes, src, sink, arcs));

        long long total = 0;
        vector<long long> net(nodes, 0);
        vector<array<long long, 3>> residual;
        for (auto& a : arcs){
            long long f = g->flow[a.id];
            assert(0 <= f && f <= a.cap && g->flow[a.id ^ 1] == -f);
            net[a.u] -= f, net[a.v] += f, total += f * a.cost;
            if (f < a.cap) residual.push_back({a.u, a.v, a.cost});
            if (f > 0) residual.push_back({a.v, a.u, -a.cost});
        }
        for (int x = 0; x < nodes; x++) assert(net[x] == (x == sink ? flow : x == src ? -flow : 0));
        assert(cost == total);
        assert(!has_negative_cycle(nodes, residual));
    }
    delete g;
    return 0;
}
