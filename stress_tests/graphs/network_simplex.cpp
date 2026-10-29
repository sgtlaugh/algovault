#include "../common.h"

#define main library_main
#include "../../code_library/graphs/network_simplex.cpp"
#undef main

#define main circulation_main
#include "../../code_library/graphs/min_cost_circulation.cpp"
#undef main

struct Arc{ int u, v; long long low, high, cost; };

struct Instance{
    int n;
    vector<long long> supply;
    vector<Arc> arcs;
};

/// Tries every integer flow vector within the bounds, returns the min cost or nullopt
optional<long long> brute_force(const Instance& g){
    int m = g.arcs.size();
    vector<long long> f(m);
    for (int i = 0; i < m; i++) f[i] = g.arcs[i].low;

    optional<long long> best;
    while (true){
        vector<long long> net = g.supply;
        long long cost = 0;
        for (int i = 0; i < m; i++) net[g.arcs[i].u] -= f[i], net[g.arcs[i].v] += f[i], cost += f[i] * g.arcs[i].cost;
        bool balanced = true;
        for (long long x : net) balanced &= x == 0;
        if (balanced && (!best || cost < *best)) best = cost;

        int i = 0;
        while (i < m && f[i] == g.arcs[i].high) f[i] = g.arcs[i].low, i++;
        if (i == m) return best;
        f[i]++;
    }
}

/// The b-flow as a circulation through a source S and sink T, every unit of supply paid -big to force it through
optional<long long> by_circulation(const Instance& g, bool zero_costs){
    int n = g.n, src = n, snk = n + 1;
    MinCostCirculation c(n + 2);
    vector<long long> b = g.supply;
    long long big = 1, shift = 0, total = 0;
    for (auto& a : g.arcs){
        long long cost = zero_costs ? 0 : a.cost;
        c.add_edge(a.u, a.v, a.high - a.low, cost);
        b[a.u] -= a.low, b[a.v] += a.low, shift += a.low * cost, big += abs(cost);
    }

    long long balance = 0;
    for (int v = 0; v < n; v++){
        balance += b[v];
        if (b[v] > 0) c.add_edge(src, v, b[v], 0), total += b[v];
        if (b[v] < 0) c.add_edge(v, snk, -b[v], 0);
    }
    if (balance != 0) return nullopt;

    int back = c.add_edge(snk, src, total, -big);
    long long cost = c.solve();
    if (c.flow(back) < total) return nullopt;
    return cost + big * total + shift;
}

/// Checks bounds, conservation, the returned cost, complementary slackness (an LP optimality proof) and the dual range
void check_certificate(const Instance& g, NetworkSimplex& ns, const vector<int>& ids, __int128 cost){
    long long max_cost = 1;
    vector<long long> net = g.supply;
    __int128 total = 0;
    assert((int)ns.potential.size() == g.n);
    for (int i = 0; i < (int)g.arcs.size(); i++){
        auto [u, v, low, high, c] = g.arcs[i];
        long long f = ns.flow(ids[i]);
        assert(low <= f && f <= high);
        net[u] -= f, net[v] += f, total += (__int128)f * c;
        __int128 rc = (__int128)c + ns.potential[u] - ns.potential[v];
        if (f < high) assert(rc >= 0);
        if (f > low) assert(rc <= 0);
        max_cost = max(max_cost, abs(c));
    }

    for (long long x : net) assert(x == 0);
    assert(cost == total);
    for (long long p : ns.potential) assert(abs(p) <= 2 * g.n * max_cost + 1);
}

Instance random_instance(int n, int m, long long max_bound, long long max_width, long long max_cost, long long max_supply){
    Instance g{n, vector<long long>(n), {}};
    for (int i = 0; i < m; i++){
        long long low = stress::rand_int(-max_bound, max_bound), high = min(max_bound, low + stress::rand_int(0, max_width));
        g.arcs.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), low, high, stress::rand_int(-max_cost, max_cost)});
    }

    /// Mostly balanced supplies, sometimes off by one so the total supply check is exercised
    long long sum = 0;
    for (int v = 0; v + 1 < n; v++) g.supply[v] = stress::rand_int(-max_supply, max_supply), sum += g.supply[v];
    g.supply[n - 1] = -sum + (stress::rand_int(0, 9) == 0 ? stress::rand_int(-1, 1) : 0);
    return g;
}

/// Supplies read off a random flow within the bounds, so the instance is feasible however wide the values are
Instance feasible(Instance g){
    fill(g.supply.begin(), g.supply.end(), 0);
    for (auto& a : g.arcs){
        long long f = stress::rand_int(a.low, a.high);
        g.supply[a.u] += f, g.supply[a.v] -= f;
    }

    return g;
}

pair<optional<__int128>, vector<int>> run(const Instance& g, NetworkSimplex& ns){
    vector<int> ids;
    for (int v = 0; v < g.n; v++){
        if (g.supply[v] != 0) ns.add_supply(v, g.supply[v]);
    }
    for (auto& a : g.arcs){
        if (a.low == 0 && stress::rand_int(0, 1)) ids.push_back(ns.add_edge(a.u, a.v, a.high, a.cost));
        else ids.push_back(ns.add_edge(a.u, a.v, a.low, a.high, a.cost));
    }
    auto cost = ns.solve();
    assert((int)ns.potential.size() == g.n);

    /// Strong feasibility, which rules out cycling: every node can push flow up to the root, a wrong tie rule breaks it
    for (int x = 0; x < g.n; x++) assert(ns.residual(x, false) > 0);
    return {cost, ids};
}

void check_large(const Instance& g){
    NetworkSimplex ns(g.n);
    auto [cost, ids] = run(g, ns);
    assert(cost.has_value() == by_circulation(g, true).has_value());
    if (cost) check_certificate(g, ns, ids, *cost);
}

/// Brute force on tiny graphs, the circulation and the certificate on medium ones, judge limits and degenerate ties
int main(){
    NetworkSimplex empty(0);
    assert(empty.solve() == 0);

    /// 1e12 units on one edge of cost 1e9: a single flow * cost product already overflows long long
    NetworkSimplex heavy(2);
    heavy.add_supply(0, 1000000000000);
    heavy.add_supply(1, -1000000000000);
    heavy.add_edge(0, 1, 2000000000000, 1000000000);
    assert(heavy.solve() == (__int128)1000000000000 * 1000000000);

    /// Negative cost edges start full, moving 4e18 units out of node 0 at once, at the header's sum (high - low) limit
    Instance star{5, vector<long long>(5), {}};
    for (int i = 1; i < 5; i++){
        star.arcs.push_back({0, i, 0, 1000000000000000000, -1});
        star.arcs.push_back({i, 0, 0, 1000000000000000000, 0});
    }
    NetworkSimplex saturated(5);
    auto [star_cost, star_ids] = run(star, saturated);
    assert(star_cost == -4000000000000000000);
    check_certificate(star, saturated, star_ids, *star_cost);

    for (long long it = 0; it < stress::scaled(2500); it++){
        int n = stress::rand_int(1, 4), m = stress::rand_int(0, it % 2 ? 4 : 6);
        Instance g = random_instance(n, m, 2, it % 2 ? 4 : 2, 5, 3);
        NetworkSimplex ns(n);
        auto [cost, ids] = run(g, ns);
        auto best = brute_force(g);
        assert(cost.has_value() == best.has_value());
        if (!cost) continue;
        assert(*cost == *best);
        check_certificate(g, ns, ids, *cost);
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(1, 30), m = stress::rand_int(0, 5 * n);
        /// Degenerate classes: unit caps with costs in {-1, 0, 1} make many ties for the leaving arc rule
        bool ties = it % 3 == 0;
        Instance g = ties ? random_instance(n, m, 1, 1, 1, 2) : random_instance(n, m, 1000, 1000, 10000, 3000);
        if (it % 2) g = feasible(g);
        NetworkSimplex ns(n);
        auto [cost, ids] = run(g, ns);
        auto expected = by_circulation(g, false);
        assert(cost.has_value() == expected.has_value());
        if (!cost) continue;
        assert(*cost == *expected);
        check_certificate(g, ns, ids, *cost);
    }

    /// Edges added after a solve: the second solve must equal a fresh one on the full graph
    for (long long it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(1, 12);
        Instance g = random_instance(n, stress::rand_int(0, 4 * n), 50, 50, 100, 50);
        Instance extra = random_instance(n, stress::rand_int(1, 2 * n), 50, 50, 100, 0);
        NetworkSimplex ns(n);
        auto [first, ids] = run(g, ns);
        for (auto& a : extra.arcs) ids.push_back(ns.add_edge(a.u, a.v, a.low, a.high, a.cost)), g.arcs.push_back(a);
        auto again = ns.solve();
        NetworkSimplex fresh(n);
        assert(again == run(g, fresh).first);
        if (again) check_certificate(g, ns, ids, *again);
    }

    /// Library Checker limits: n = 100, m = 1000, every |value| up to 1e9, so the total cost overflows long long
    for (long long it = 0; it < stress::scaled(20); it++){
        long long width = it % 2 ? 2000000000 : 1000;
        Instance g = random_instance(100, 1000, 1000000000, width, 1000000000, 1000000000);
        check_large(it % 4 < 2 ? feasible(g) : g);
    }

    /// n max|cost| = 1e18, the header's limit, where reduced costs near long long's range
    for (long long it = 0; it < stress::scaled(200); it++){
        int n = 10;
        Instance g = random_instance(n, stress::rand_int(0, 40), 100, 100, 100000000000000000, 300);
        check_large(it % 2 ? feasible(g) : g);
    }

    /// A long negative ring at |cost| = 1e9 and supplies that need long paths, so the tree gets deep
    const int n = 2000;
    Instance g{n, vector<long long>(n), {}};
    for (int i = 0; i < n; i++) g.arcs.push_back({i, (i + 1) % n, 0, stress::rand_int(1, 1000000000), -1000000000});
    for (int i = 0; i < 4 * n; i++){
        g.arcs.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), 0, stress::rand_int(0, 1000000000), stress::rand_int(-1000000000, 1000000000)});
    }
    shuffle(g.arcs.begin(), g.arcs.end(), stress::rng());
    check_large(feasible(g));

    return 0;
}
