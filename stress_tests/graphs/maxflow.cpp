#include "../common.h"

#define main library_main
#include "../../code_library/graphs/maxflow.cpp"
#undef main

struct Arc{ int u, v; long long cap; };

long long cut_value(int mask, const vector<Arc>& arcs){
    long long cut = 0;
    for (auto& a : arcs) if ((mask >> a.u & 1) && !(mask >> a.v & 1)) cut += a.cap;
    return cut;
}

/// Max flow = min cut, enumerated over every source side containing src but not sink
long long min_cut(int n, int src, int sink, const vector<Arc>& arcs){
    long long best = LLONG_MAX;
    for (int mask = 0; mask < (1 << n); mask++){
        if (!(mask >> src & 1) || (mask >> sink & 1)) continue;
        best = min(best, cut_value(mask, arcs));
    }
    return best;
}

/// The residual-reachable side is the intersection of every min cut source side
int smallest_min_cut_side(int n, int src, int sink, const vector<Arc>& arcs){
    long long best = min_cut(n, src, sink, arcs);
    int side = (1 << n) - 1;
    for (int mask = 0; mask < (1 << n); mask++){
        if ((mask >> src & 1) && !(mask >> sink & 1) && cut_value(mask, arcs) == best) side &= mask;
    }
    return side;
}

int source_side_mask(int n, const function<bool(int)>& in_source_side){
    int side = 0;
    for (int x = 0; x < n; x++) if (in_source_side(x)) side |= 1 << x;
    return side;
}

/// Source side against the enumerated min cuts, cut edges against the arcs crossing it
void check_cut(const FlowGraph& g, long long value, const vector<Arc>& arcs){
    int side = source_side_mask(g.n, [&](int x){ return g.in_source_side(x); });
    assert(side == smallest_min_cut_side(g.n, g.src, g.sink, arcs));

    long long total = 0;
    auto ids = g.cut_edges();
    for (int id : ids){
        auto& e = g.E[id];
        assert(e.cap > 0 && (side >> e.u & 1) && !(side >> e.v & 1));
        total += e.cap;
    }

    long long crossing = count_if(arcs.begin(), arcs.end(), [&](const Arc& a){ return a.cap > 0 && (side >> a.u & 1) && !(side >> a.v & 1); });
    assert(total == value && (long long)ids.size() == crossing);
}

/// Capacity respected on every edge, conservation everywhere but src and sink
void check_flow(const FlowGraph& g, long long value){
    vector<long long> net(g.n, 0);
    for (size_t id = 0; id < g.E.size(); id += 2){
        auto& e = g.E[id];
        assert(0 <= e.flow && e.flow <= e.cap && g.E[id + 1].flow == -e.flow);
        net[e.u] -= e.flow, net[e.v] += e.flow;
    }
    for (int x = 0; x < g.n; x++){
        if (x == g.src) assert(net[x] == -value);
        else if (x == g.sink) assert(net[x] == value);
        else assert(net[x] == 0);
    }
}

/// Header recipe: max closure by min cut against every subset that is closed under the forcing pairs
void check_closure(const vector<long long>& weight, const vector<pair<int, int>>& forces){
    int n = weight.size(), src = n, sink = n + 1;
    FlowGraph g(n + 2, src, sink);
    long long positive = 0;
    for (int i = 0; i < n; i++){
        if (weight[i] > 0) g.add_directed_edge(src, i, weight[i]), positive += weight[i];
        else if (weight[i] < 0) g.add_directed_edge(i, sink, -weight[i]);
    }
    for (auto [i, j] : forces) g.add_directed_edge(i, j, LLONG_MAX);
    long long value = positive - g.maxflow();

    auto closed = [&](int mask){
        for (auto [i, j] : forces) if ((mask >> i & 1) && !(mask >> j & 1)) return false;
        return true;
    };
    auto total = [&](int mask){
        long long sum = 0;
        for (int i = 0; i < n; i++) if (mask >> i & 1) sum += weight[i];
        return sum;
    };
    long long best = 0;
    for (int mask = 0; mask < (1 << n); mask++) if (closed(mask)) best = max(best, total(mask));

    int picked = source_side_mask(n, [&](int x){ return g.in_source_side(x); });
    assert(value == best && closed(picked) && total(picked) == best);
}

/// Header recipe: Dinkelbach over closures, returns the densest vertex set it settles on
vector<int> densest_subgraph(int n, const vector<pair<int, int>>& edges){
    int m = edges.size(), src = n + m, sink = n + m + 1;
    vector<int> chosen(n);
    iota(chosen.begin(), chosen.end(), 0);
    long long a = m, b = n;

    while (true){
        FlowGraph g(n + m + 2, src, sink);
        for (int i = 0; i < m; i++){
            g.add_directed_edge(src, n + i, b);
            g.add_directed_edge(n + i, edges[i].first, LLONG_MAX), g.add_directed_edge(n + i, edges[i].second, LLONG_MAX);
        }
        for (int v = 0; v < n; v++) g.add_directed_edge(v, sink, a);
        if (m * b - g.maxflow() <= 0) return chosen;

        chosen.clear();
        for (int v = 0; v < n; v++) if (g.in_source_side(v)) chosen.push_back(v);
        a = count_if(edges.begin(), edges.end(), [&](pair<int, int> e){ return g.in_source_side(e.first) && g.in_source_side(e.second); });
        b = chosen.size();
    }
}

/// Densest subgraph against the best |E(S)| / |S| over every nonempty vertex set, compared as fractions
void check_density(int n, const vector<pair<int, int>>& edges){
    auto induced = [&](int mask){
        return count_if(edges.begin(), edges.end(), [&](pair<int, int> e){ return (mask >> e.first & 1) && (mask >> e.second & 1); });
    };
    long long best_e = 0, best_v = 1;
    for (int mask = 1; mask < (1 << n); mask++){
        long long e = induced(mask), v = __builtin_popcount(mask);
        if (e * best_v > best_e * v) best_e = e, best_v = v;
    }

    int mask = 0;
    for (int v : densest_subgraph(n, edges)) mask |= 1 << v;
    assert(mask && induced(mask) * best_v == best_e * __builtin_popcount(mask));
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(2, 8), src = stress::rand_int(0, n - 1), sink = (src + stress::rand_int(1, n - 1)) % n;
        long long max_cap = it % 3 ? 10 : 1000000000000000LL;

        FlowGraph g(n, src, sink);
        DenseFlowGraph dense(n, src, sink);
        vector<Arc> arcs;
        for (int m = stress::rand_int(0, 20); m; m--){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long cap = stress::rand_int(0, max_cap);
            if (stress::rand_int(0, 1)){
                g.add_directed_edge(u, v, cap), dense.add_directed_edge(u, v, cap);
                arcs.push_back({u, v, cap});
            }
            else{
                g.add_edge(u, v, cap), dense.add_edge(u, v, cap);
                arcs.push_back({u, v, cap}), arcs.push_back({v, u, cap});
            }
        }

        assert(g.cut_edges().empty());
        assert(source_side_mask(n, [&](int x){ return g.in_source_side(x) || dense.in_source_side(x); }) == 0);

        long long flow = g.maxflow();
        assert(flow == min_cut(n, src, sink, arcs));
        check_flow(g, flow);
        check_cut(g, flow, arcs);
        assert(dense.maxflow() == flow);
        assert(source_side_mask(n, [&](int x){ return dense.in_source_side(x); }) == source_side_mask(n, [&](int x){ return g.in_source_side(x); }));
    }

    /// Found by search: reaching flow 3 here needs a reverse edge to cancel earlier flow, without one both get 2
    {
        vector<Arc> arcs = {{0, 2, 1}, {0, 1, 1}, {1, 3, 3}, {0, 3, 1}, {3, 2, 2}, {3, 5, 2}, {4, 5, 1}, {4, 2, 1}, {1, 4, 2}, {3, 0, 3}, {2, 3, 2}};
        FlowGraph g(6, 0, 5);
        DenseFlowGraph dense(6, 0, 5);
        for (auto& a : arcs) g.add_directed_edge(a.u, a.v, a.cap), dense.add_directed_edge(a.u, a.v, a.cap);
        assert(min_cut(6, 0, 5, arcs) == 3 && g.maxflow() == 3 && dense.maxflow() == 3);
    }

    /// Dense graphs: the matrix version against the edge-list version
    for (long long it = 0; it < stress::scaled(4); it++){
        int n = 400;
        FlowGraph g(n, 0, n - 1);
        DenseFlowGraph dense(n, 0, n - 1);
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++){
                if (u == v || stress::rand_int(0, 1)) continue;
                long long cap = stress::rand_int(1, 1000000000);
                g.add_directed_edge(u, v, cap), dense.add_directed_edge(u, v, cap);
            }
        }
        assert(dense.maxflow() == g.maxflow());
        for (int x = 0; x < n; x++) assert(dense.in_source_side(x) == g.in_source_side(x));
    }

    /// Node capacities limit every node, src and sink included: split x into x_in = 2x and x_out = 2x + 1
    for (long long it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(2, 6), src = stress::rand_int(0, n - 1), sink = (src + stress::rand_int(1, n - 1)) % n;
        vector<long long> node_cap(n);
        for (auto& c : node_cap) c = stress::rand_int(0, 12);

        vector<Arc> arcs;
        for (int x = 0; x < n; x++) arcs.push_back({2 * x, 2 * x + 1, node_cap[x]});
        FlowGraphWithNodeCap g(n, src, sink, node_cap);
        for (int m = stress::rand_int(0, 14); m; m--){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            long long cap = stress::rand_int(0, 10);
            if (stress::rand_int(0, 1)){
                g.add_directed_edge(u, v, cap);
                arcs.push_back({2 * u + 1, 2 * v, cap});
            }
            else{
                g.add_edge(u, v, cap);
                arcs.push_back({2 * u + 1, 2 * v, cap}), arcs.push_back({2 * v + 1, 2 * u, cap});
            }
        }

        long long flow = g.maxflow();
        assert(flow == min_cut(2 * n, 2 * src, 2 * sink + 1, arcs));
        check_flow(g.flowgraph, flow);
        check_cut(g.flowgraph, flow, arcs);
    }

    /// Node splitting doubles the node count, past the old fixed MAXN = 50010 arrays at n > 25005
    {
        int n = 30000;
        vector<long long> node_cap(n, 1);
        node_cap[0] = node_cap[1] = n;
        FlowGraphWithNodeCap g(n, 0, 1, node_cap);
        for (int v = 2; v < n; v++) g.add_directed_edge(0, v, 1), g.add_directed_edge(v, 1, 1);
        assert(g.maxflow() == n - 2);
    }

    /// Max closure: small weights, then weights near 1e12 so LLONG_MAX forcing edges sit beside large finite ones
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 10);
        long long max_weight = it % 4 ? 10 : 1000000000000LL;
        vector<long long> weight(n);
        for (auto& w : weight) w = stress::rand_int(-max_weight, max_weight);
        vector<pair<int, int>> forces;
        for (int m = stress::rand_int(0, 15); m; m--) forces.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1)});
        check_closure(weight, forces);
    }

    /// Max density subgraph on simple graphs, edgeless ones included
    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, 8);
        vector<pair<int, int>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) if (stress::rand_int(0, 2) == 0) edges.push_back({u, v});
        }
        check_density(n, edges);
    }

    /// K4 with a pendant vertex: the whole graph has density 7 / 5, the K4 alone 6 / 4
    assert((densest_subgraph(5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}, {0, 4}}) == vector<int>{0, 1, 2, 3}));

    /// Unit-capacity bipartite matching with a planted perfect matching, the O(E sqrt(V)) case
    {
        int half = 20000, src = 2 * half, sink = 2 * half + 1;
        FlowGraph g(2 * half + 2, src, sink);
        for (int i = 0; i < half; i++){
            g.add_directed_edge(src, i, 1), g.add_directed_edge(half + i, sink, 1);
            for (int k = 0; k < 5; k++) g.add_directed_edge(i, half + stress::rand_int(0, half - 1), 1);
            g.add_directed_edge(i, half + i, 1);
        }
        assert(g.maxflow() == half);
    }
    return 0;
}
