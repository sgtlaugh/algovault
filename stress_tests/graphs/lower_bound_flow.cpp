#include "../common.h"

#define main library_main
#include "../../code_library/graphs/lower_bound_flow.cpp"
#undef main

const long long NONE = LowerBoundFlow::INFEASIBLE;

/// Every edge within its bounds, conservation everywhere but s and t, and the net flow out of s equals value
void check_flow(LowerBoundFlow& g, const vector<LowerBoundFlow::Bound>& edges, int s, int t, long long value){
    vector<long long> net(g.n, 0);
    for (int id = 0; id < (int)edges.size(); id++){
        long long f = g.flow(id);
        assert(edges[id].low <= f && f <= edges[id].high);
        net[edges[id].u] += f, net[edges[id].v] -= f;
    }
    for (int x = 0; x < g.n; x++){
        if (x == s) assert(net[x] == value);
        else if (x == t) assert(net[x] == -value);
        else assert(net[x] == 0);
    }
}

LowerBoundFlow make_graph(int n, const vector<LowerBoundFlow::Bound>& edges){
    LowerBoundFlow g(n);
    for (int id = 0; id < (int)edges.size(); id++) assert(g.add_edge(edges[id].u, edges[id].v, edges[id].low, edges[id].high) == id);
    return g;
}

/// Exhaustive over every integer flow: per (s, t) the min and max value, NONE when no flow conserves at the other nodes
struct Brute{
    bool circulation = false;
    vector<vector<long long>> lo, hi;

    Brute(int n, const vector<LowerBoundFlow::Bound>& edges) : lo(n, vector<long long>(n, LLONG_MAX)), hi(n, vector<long long>(n, LLONG_MIN)){
        int m = edges.size();
        vector<long long> f(m);
        for (int i = 0; i < m; i++) f[i] = edges[i].low;

        while (true){
            vector<long long> net(n, 0);
            for (int i = 0; i < m; i++) net[edges[i].u] += f[i], net[edges[i].v] -= f[i];

            vector<int> bad;
            for (int x = 0; x < n; x++) if (net[x]) bad.push_back(x);
            if (bad.empty()) circulation = true;
            for (int s = 0; s < n; s++){
                for (int t = 0; t < n; t++){
                    if (s == t) continue;
                    bool ok = true;
                    for (int x : bad) ok &= x == s || x == t;
                    if (ok) lo[s][t] = min(lo[s][t], net[s]), hi[s][t] = max(hi[s][t], net[s]);
                }
            }

            int i = 0;
            while (i < m && f[i] == edges[i].high) f[i] = edges[i].low, i++;
            if (i == m) break;
            f[i]++;
        }
    }
};

/// Hoffman: a circulation exists iff no node set S takes in more forced flow than it can send out
bool cut_feasible(int n, const vector<LowerBoundFlow::Bound>& edges, int s, int t){
    for (int mask = 0; mask < (1 << n); mask++){
        if (s != -1 && (mask >> s & 1) != (mask >> t & 1)) continue;
        long long in_low = 0, out_high = 0;
        for (auto& e : edges){
            bool from = mask >> e.u & 1, to = mask >> e.v & 1;
            if (!from && to) in_low += e.low;
            if (from && !to) out_high += e.high;
        }
        if (in_low > out_high) return false;
    }
    return true;
}

/// Max flow = min over s-t cuts of high out minus low in, min flow = max over s-t cuts of low out minus high in
pair<long long, long long> cut_bounds(int n, const vector<LowerBoundFlow::Bound>& edges, int s, int t){
    long long mn = LLONG_MIN, mx = LLONG_MAX;
    for (int mask = 0; mask < (1 << n); mask++){
        if (!(mask >> s & 1) || (mask >> t & 1)) continue;
        long long out_high = 0, out_low = 0, in_high = 0, in_low = 0;
        for (auto& e : edges){
            bool from = mask >> e.u & 1, to = mask >> e.v & 1;
            if (from && !to) out_high += e.high, out_low += e.low;
            if (!from && to) in_high += e.high, in_low += e.low;
        }
        mx = min(mx, out_high - in_low), mn = max(mn, out_low - in_high);
    }
    return {mn, mx};
}

vector<LowerBoundFlow::Bound> random_edges(int n, int m, long long max_cap){
    vector<LowerBoundFlow::Bound> edges;
    for (int i = 0; i < m; i++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        long long x = stress::rand_int(0, max_cap), y = stress::rand_int(0, max_cap);
        if (stress::rand_int(0, 2) == 0) x = 0;
        edges.push_back({u, v, min(x, y), max(x, y)});
    }
    return edges;
}

/// Random circulation as a sum of random cycles, then bounds around it so a feasible flow is planted
vector<LowerBoundFlow::Bound> planted_edges(int n, int cycles, long long max_cap){
    vector<LowerBoundFlow::Bound> edges;
    vector<long long> f;
    for (int c = 0; c < cycles; c++){
        int len = stress::rand_int(2, min(n, 8));
        long long amount = stress::rand_int(0, max_cap);
        vector<int> nodes;
        for (int i = 0; i < len; i++) nodes.push_back(stress::rand_int(0, n - 1));
        for (int i = 0; i < len; i++) edges.push_back({nodes[i], nodes[(i + 1) % len], 0, 0}), f.push_back(amount);
    }
    for (size_t i = 0; i < edges.size(); i++){
        edges[i].low = stress::rand_int(0, f[i]), edges[i].high = f[i] + stress::rand_int(0, max_cap);
    }
    return edges;
}

void check_all(int n, const vector<LowerBoundFlow::Bound>& edges, int s, int t, long long want_min, long long want_max){
    LowerBoundFlow g = make_graph(n, edges);

    long long any = g.feasible_flow(s, t);
    if (want_max == NONE){
        assert(any == NONE && g.max_flow(s, t) == NONE && g.min_flow(s, t) == NONE);
        return;
    }
    assert(want_min <= any && any <= want_max);
    check_flow(g, edges, s, t, any);

    assert(g.max_flow(s, t) == want_max);
    check_flow(g, edges, s, t, want_max);
    assert(g.min_flow(s, t) == want_min);
    check_flow(g, edges, s, t, want_min);
}

int main(){
    for (long long it = 0; it < stress::scaled(2500); it++){
        int n = stress::rand_int(2, 4), m = stress::rand_int(0, 5);
        auto edges = random_edges(n, m, 3);
        Brute brute(n, edges);

        LowerBoundFlow g = make_graph(n, edges);
        assert(g.circulation() == brute.circulation);
        if (brute.circulation) check_flow(g, edges, 0, 1, 0);

        for (int s = 0; s < n; s++){
            for (int t = 0; t < n; t++){
                if (s == t) continue;
                bool feasible = brute.hi[s][t] != LLONG_MIN;
                check_all(n, edges, s, t, feasible ? brute.lo[s][t] : NONE, feasible ? brute.hi[s][t] : NONE);
            }
        }
    }

    /// Larger capacities on up to 10 nodes, against Hoffman's condition and the cut formulas
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(2, 10), s = stress::rand_int(0, n - 1), t = (s + stress::rand_int(1, n - 1)) % n;
        long long max_cap = it % 3 ? 20 : 1000000000000LL;
        auto edges = it % 2 ? random_edges(n, stress::rand_int(0, 25), max_cap) : planted_edges(n, stress::rand_int(1, 4), max_cap);

        LowerBoundFlow g = make_graph(n, edges);
        bool circ = cut_feasible(n, edges, -1, -1);
        assert(g.circulation() == circ);
        if (circ) check_flow(g, edges, s, t, 0);

        if (!cut_feasible(n, edges, s, t)){
            check_all(n, edges, s, t, NONE, NONE);
            continue;
        }
        auto [mn, mx] = cut_bounds(n, edges, s, t);
        check_all(n, edges, s, t, mn, mx);
    }

    /// Big planted instances: certificates only, every answer is a valid flow and planted circulations stay feasible
    for (long long it = 0; it < stress::scaled(6); it++){
        int n = 3000;
        auto edges = planted_edges(n, 4000, 1000000000);
        for (int i = 0; i < 3000; i++){
            long long high = stress::rand_int(0, 1000000000);
            edges.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), 0, high});
        }
        LowerBoundFlow g = make_graph(n, edges);
        assert(g.circulation());
        check_flow(g, edges, 0, 1, 0);

        int s = stress::rand_int(0, n - 1), t = (s + stress::rand_int(1, n - 1)) % n;
        long long mx = g.max_flow(s, t);
        check_flow(g, edges, s, t, mx);
        long long mn = g.min_flow(s, t);
        check_flow(g, edges, s, t, mn);
        assert(mn <= 0 && 0 <= mx);
    }

    /// A path of 100000 nodes is a level graph 100000 deep; one shared lower bound keeps the S -> T phase short, random ones make Dinic quadratic here
    {
        int n = 100000;
        vector<LowerBoundFlow::Bound> edges;
        long long low = stress::rand_int(0, 1000), min_high = LLONG_MAX;
        for (int x = 0; x + 1 < n; x++){
            long long high = stress::rand_int(1000, 1000000000);
            edges.push_back({x, x + 1, low, high});
            min_high = min(min_high, high);
        }
        check_all(n, edges, 0, n - 1, low, min_high);
        edges[n / 2].low = edges[n / 2].high = 2000000000;
        check_all(n, edges, 0, n - 1, NONE, NONE);
    }

    /// Upper bounds summing to exactly LLONG_MAX, the documented limit, with every internal sum at the edge of overflow
    {
        vector<LowerBoundFlow::Bound> pinned = {{1, 0, LLONG_MAX, LLONG_MAX}};
        check_all(2, pinned, 0, 1, -LLONG_MAX, -LLONG_MAX);
        check_all(2, pinned, 1, 0, LLONG_MAX, LLONG_MAX);

        long long third = LLONG_MAX / 3;
        vector<LowerBoundFlow::Bound> cycle = {{0, 1, 1, third}, {1, 2, 2, third}, {2, 0, 3, LLONG_MAX - 2 * third}};
        assert(make_graph(3, cycle).circulation());
        for (int s = 0; s < 3; s++){
            for (int t = 0; t < 3; t++){
                if (s == t) continue;
                auto [mn, mx] = cut_bounds(3, cycle, s, t);
                check_all(3, cycle, s, t, mn, mx);
            }
        }
    }
    return 0;
}
