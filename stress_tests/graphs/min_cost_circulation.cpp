#include "../common.h"

#define main library_main
#include "../../code_library/graphs/min_cost_circulation.cpp"
#undef main

struct Arc{ int u, v; long long cap, cost; };

/// Tries every integer flow vector, returns {max flow from s to t, min cost at that flow}, s = t = -1 asks for a circulation
pair<long long, long long> brute_force(int n, const vector<Arc>& arcs, int s, int t){
    int m = arcs.size();
    vector<long long> f(m, 0);
    pair<long long, long long> best = {-1, 0};
    while (true){
        vector<long long> net(n, 0);
        long long cost = 0;
        for (int i = 0; i < m; i++) net[arcs[i].u] -= f[i], net[arcs[i].v] += f[i], cost += f[i] * arcs[i].cost;
        bool balanced = true;
        for (int x = 0; x < n; x++) balanced &= x == s || x == t || net[x] == 0;
        long long value = s == -1 ? 0 : net[t];
        if (balanced && (value > best.first || (value == best.first && cost < best.second))) best = {value, cost};

        int i = 0;
        while (i < m && f[i] == arcs[i].cap) f[i++] = 0;
        if (i == m) return best;
        f[i]++;
    }
}

/// A circulation is optimal exactly when its residual graph has no negative cycle
bool has_negative_cycle(int n, const vector<array<long long, 3>>& residual){
    vector<long long> dist(n, 0);
    for (int round = 0; round <= n; round++){
        bool relaxed = false;
        for (auto& [u, v, w] : residual){
            if (dist[u] + w < dist[v]) dist[v] = dist[u] + w, relaxed = true;
        }
        if (!relaxed) return false;
    }
    return true;
}

/// Checks capacities, conservation, the returned cost and the absence of negative residual cycles, returns the flows
vector<long long> check_optimal(int n, const vector<Arc>& arcs, MinCostCirculation& g, const vector<int>& ids, long long cost){
    long long total = 0;
    vector<long long> net(n, 0), flows;
    vector<array<long long, 3>> residual;
    for (int i = 0; i < (int)arcs.size(); i++){
        auto [u, v, cap, c] = arcs[i];
        long long f = g.flow(ids[i]);
        assert(0 <= f && f <= cap);
        net[u] -= f, net[v] += f, total += f * c;
        if (f < cap) residual.push_back({u, v, c});
        if (f > 0) residual.push_back({v, u, -c});
        flows.push_back(f);
    }

    for (int x = 0; x < n; x++) assert(net[x] == 0);
    assert(cost == total);
    assert(!has_negative_cycle(n, residual));
    return flows;
}

/// Brute force on tiny graphs, the residual negative cycle certificate on larger ones and at the n = 1e4, |cost| = 1e9 limit
int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 4), m = stress::rand_int(0, 6);
        bool flow_mode = n >= 2 && it % 3 == 0;
        int s = flow_mode ? stress::rand_int(0, n - 1) : -1, t = flow_mode ? (s + stress::rand_int(1, n - 1)) % n : -1;

        vector<Arc> arcs;
        for (int i = 0; i < m; i++){
            arcs.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(0, 3), stress::rand_int(-6, 6)});
        }

        MinCostCirculation g(n);
        vector<int> ids;
        long long big = 1, out = 0;
        for (auto& a : arcs) ids.push_back(g.add_edge(a.u, a.v, a.cap, a.cost)), big += abs(a.cost), out += a.u == s ? a.cap : 0;
        int back = flow_mode ? g.add_edge(t, s, out, -big) : -1;
        long long cost = g.solve();

        auto [best_flow, best_cost] = brute_force(n, arcs, s, t);
        if (!flow_mode){
            assert(cost == best_cost);
            check_optimal(n, arcs, g, ids, cost);
            continue;
        }

        long long flow = g.flow(back);
        assert(flow == best_flow && cost + big * flow == best_cost);
        vector<long long> net(n, 0);
        for (int i = 0; i < m; i++) net[arcs[i].u] -= g.flow(ids[i]), net[arcs[i].v] += g.flow(ids[i]);
        for (int x = 0; x < n; x++) assert(net[x] == (x == t ? flow : x == s ? -flow : 0));
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 10 ? 12 : 60), m = stress::rand_int(0, 4 * n);
        /// Either large capacities or large costs, so the total cost still fits in a long long
        long long max_cap = it % 4 == 0 ? 1000000000 : it % 4 == 1 ? 1000000 : 5, max_cost = it % 4 == 0 ? 20 : 1000000000;

        vector<Arc> arcs;
        for (int i = 0; i < m; i++){
            arcs.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(0, max_cap), stress::rand_int(-max_cost, max_cost)});
        }

        MinCostCirculation g(n);
        vector<int> ids;
        for (auto& a : arcs) ids.push_back(g.add_edge(a.u, a.v, a.cap, a.cost));
        long long cost = g.solve();
        check_optimal(n, arcs, g, ids, cost);
    }

    /// A Hamiltonian ring of -1e9 arcs plus shuffled random chords: the longest negative cycle at the n and cost limits, caps <= 1e5 keep the total in range
    const int n = 10000, m = 30000;
    const long long max_cost = 1000000000;

    vector<Arc> arcs;
    for (int i = 0; i < n; i++) arcs.push_back({i, (i + 1) % n, stress::rand_int(1, 100000), -max_cost});
    while ((int)arcs.size() < m){
        arcs.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(0, 100000), stress::rand_int(-max_cost, max_cost)});
    }
    shuffle(arcs.begin(), arcs.end(), stress::rng());

    MinCostCirculation g(n);
    vector<int> ids;
    for (auto& a : arcs) ids.push_back(g.add_edge(a.u, a.v, a.cap, a.cost));
    long long cost = g.solve();
    check_optimal(n, arcs, g, ids, cost);

    return 0;
}
