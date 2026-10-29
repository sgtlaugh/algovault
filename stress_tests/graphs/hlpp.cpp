#include "../common.h"

#define main dinic_main
#include "../../code_library/graphs/maxflow.cpp"
#undef main

#define main library_main
#include "../../code_library/graphs/hlpp.cpp"
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

/// Capacity respected on every edge, conservation everywhere but src and sink, so the preflow was fully returned
template <typename T>
void check_flow(const HLPP<T>& g, long long value){
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

/// Source side equal to Dinic's, cut edges exactly the edges crossing it, summing to the flow
template <typename T>
void check_cut(const HLPP<T>& g, const FlowGraph& dinic, long long value){
    for (int x = 0; x < g.n; x++) assert(g.in_source_side(x) == dinic.in_source_side(x));
    assert(g.cut_edges() == dinic.cut_edges());

    long long total = 0;
    for (int id : g.cut_edges()) total += g.E[id].cap;
    assert(total == value);
}

/// Same flat arc layout as Dinic's adjacency lists, so a counting sort that reorders arcs is caught
template <typename T>
void check_arcs(const HLPP<T>& g, const FlowGraph& dinic){
    assert(g.start == dinic.start);
    assert(g.to == dinic.to);
    assert(g.edge_id == dinic.edge_id);
    assert(g.rev == dinic.rev);
}

/// Builds the same random graph in both, returns the arcs; undirected edges become two arcs
template <typename T>
vector<Arc> random_graph(HLPP<T>& g, FlowGraph& dinic, int m, long long max_cap){
    vector<Arc> arcs;
    for (; m; m--){
        int u = stress::rand_int(0, g.n - 1), v = stress::rand_int(0, g.n - 1);
        long long cap = stress::rand_int(0, max_cap);
        if (stress::rand_int(0, 3)){
            int id = g.add_directed_edge(u, v, cap);
            assert(dinic.add_directed_edge(u, v, cap) == id && g.E[id].u == u && g.E[id].v == v && g.E[id].cap == cap);
            arcs.push_back({u, v, cap});
        }
        else{
            g.add_undirected_edge(u, v, cap), dinic.add_undirected_edge(u, v, cap);
            arcs.push_back({u, v, cap}), arcs.push_back({v, u, cap});
        }
    }

    return arcs;
}

template <typename T>
void compare(HLPP<T>& g, FlowGraph& dinic){
    long long flow = g.maxflow();
    assert(flow == dinic.maxflow());
    check_flow(g, flow);
    check_cut(g, dinic, flow);
    check_arcs(g, dinic);
}

int main(){
    /// Tiny graphs with multi-edges and self-loops against brute force, then a second call after adding edges
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(2, 8), src = stress::rand_int(0, n - 1), sink = (src + stress::rand_int(1, n - 1)) % n;
        long long max_cap = it % 3 ? 10 : 1000000000000000LL;

        HLPP g(n, src, sink);
        FlowGraph dinic(n, src, sink);
        vector<Arc> arcs = random_graph(g, dinic, stress::rand_int(0, 20), max_cap);

        assert(g.cut_edges().empty());
        long long flow = g.maxflow();
        assert(flow == min_cut(n, src, sink, arcs));
        assert(dinic.maxflow() == flow);
        check_flow(g, flow);
        check_cut(g, dinic, flow);
        check_arcs(g, dinic);
        int side = 0;
        for (int x = 0; x < n; x++) if (g.in_source_side(x)) side |= 1 << x;
        assert(side == smallest_min_cut_side(n, src, sink, arcs));

        vector<Arc> more = random_graph(g, dinic, stress::rand_int(0, 6), max_cap);
        arcs.insert(arcs.end(), more.begin(), more.end());
        long long added = g.maxflow();
        assert(added == dinic.maxflow() && flow + added == min_cut(n, src, sink, arcs));
        check_flow(g, flow + added);
        check_cut(g, dinic, flow + added);
        check_arcs(g, dinic);
    }

    /// HLPP<int> with src capacities summing to exactly INT_MAX, the most maxflow() can push at once
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(3, 40), m = stress::rand_int(0, 200), src = 0, sink = n - 1;
        HLPP<int> g(n, src, sink);
        FlowGraph dinic(n, src, sink);
        int out = stress::rand_int(1, 4);
        long long left = INT_MAX;
        for (int k = 0; k < out; k++){
            int cap = k + 1 == out ? left : stress::rand_int(0, left);
            int v = stress::rand_int(1, n - 1);
            g.add_directed_edge(src, v, cap), dinic.add_directed_edge(src, v, cap), left -= cap;
        }
        for (; m; m--){
            int u = stress::rand_int(1, n - 1), v = stress::rand_int(1, n - 1), cap = stress::rand_int(0, INT_MAX);
            g.add_directed_edge(u, v, cap), dinic.add_directed_edge(u, v, cap);
        }
        compare(g, dinic);
    }

    /// Capacities at src summing past the max of T would overflow an excess, so maxflow() must refuse them
    stress::expect_assert("graphs/hlpp.cpp", []{
        HLPP<int> g(3, 0, 2);
        g.add_directed_edge(0, 1, INT_MAX), g.add_directed_edge(0, 1, 1), g.add_directed_edge(1, 2, 1);
        g.maxflow();
    });
    stress::expect_assert("graphs/hlpp.cpp", []{ HLPP(3, 1, 1).maxflow(); });

    /// Medium graphs against Dinic: sparse, dense, unit capacities, layered, and LLONG_MAX edges inside
    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(2, 150), kind = it % 5;
        int src = stress::rand_int(0, n - 1), sink = (src + stress::rand_int(1, n - 1)) % n;
        HLPP g(n, src, sink);
        FlowGraph dinic(n, src, sink);

        if (kind == 0) random_graph(g, dinic, stress::rand_int(0, 3 * n), 1000000000);
        else if (kind == 1) random_graph(g, dinic, stress::rand_int(0, n * n), 100);
        else if (kind == 2) random_graph(g, dinic, stress::rand_int(0, 6 * n), 1);
        else if (kind == 3){
            int layers = stress::rand_int(2, 6);
            auto layer = [&](int x){ return x == src ? 0 : x == sink ? layers + 1 : 1 + x % layers; };
            for (int e = stress::rand_int(0, 8 * n); e; e--){
                int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
                if (layer(v) != layer(u) + 1 && stress::rand_int(0, 9)) continue;
                long long cap = stress::rand_int(1, 20);
                g.add_directed_edge(u, v, cap), dinic.add_directed_edge(u, v, cap);
            }
        }
        else{
            for (int e = stress::rand_int(0, 4 * n); e; e--){
                int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
                long long cap = stress::rand_int(0, 1) ? LLONG_MAX : stress::rand_int(0, 1000000000000000000LL);
                if (u == src || v == src) cap = stress::rand_int(0, LLONG_MAX / (4 * n));
                g.add_directed_edge(u, v, cap), dinic.add_directed_edge(u, v, cap);
            }
        }
        compare(g, dinic);
    }

    /// Unit-capacity bipartite matching with a planted perfect matching, many gaps and global relabels
    {
        int half = 20000, src = 2 * half, sink = 2 * half + 1;
        HLPP g(2 * half + 2, src, sink);
        FlowGraph dinic(2 * half + 2, src, sink);
        for (int i = 0; i < half; i++){
            g.add_directed_edge(src, i, 1), g.add_directed_edge(half + i, sink, 1);
            dinic.add_directed_edge(src, i, 1), dinic.add_directed_edge(half + i, sink, 1);
            for (int k = 0; k < 4; k++){
                int j = half + stress::rand_int(0, half - 1);
                g.add_directed_edge(i, j, 1), dinic.add_directed_edge(i, j, 1);
            }
            g.add_directed_edge(i, half + i, 1), dinic.add_directed_edge(i, half + i, 1);
        }
        compare(g, dinic);
    }

    /// Bottlenecked fan: most of the preflow stalls behind a thin cut and must flow back to src in the second pass
    {
        int n = 30000, src = 0, sink = n - 1;
        HLPP g(n, src, sink);
        FlowGraph dinic(n, src, sink);
        for (int v = 1; v < n - 2; v++){
            long long cap = stress::rand_int(1, 1000000);
            g.add_directed_edge(src, v, cap), dinic.add_directed_edge(src, v, cap);
            int w = stress::rand_int(1, n - 3);
            g.add_directed_edge(v, w, cap), dinic.add_directed_edge(v, w, cap);
            g.add_directed_edge(v, n - 2, 1), dinic.add_directed_edge(v, n - 2, 1);
        }
        g.add_directed_edge(n - 2, sink, 1000), dinic.add_directed_edge(n - 2, sink, 1000);
        compare(g, dinic);
    }

    /// A 200000 node path: heights reach n - 1, so the level lists and active stacks span every height
    {
        int n = 200000;
        HLPP g(n, 0, n - 1);
        long long bottleneck = LLONG_MAX;
        for (int i = 0; i + 1 < n; i++){
            long long cap = stress::rand_int(500, 1000);
            g.add_directed_edge(i, i + 1, cap), bottleneck = min(bottleneck, cap);
        }
        assert(g.maxflow() == bottleneck);
        check_flow(g, bottleneck);
    }

    return 0;
}
