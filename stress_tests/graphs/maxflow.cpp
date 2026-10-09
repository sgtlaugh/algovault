#include "../common.h"

#define main library_main
#include "../../code_library/graphs/maxflow.cpp"
#undef main

struct Arc{ int u, v; long long cap; };

/// Max flow = min cut, enumerated over every source side containing src but not sink
long long min_cut(int n, int src, int sink, const vector<Arc>& arcs){
    long long best = LLONG_MAX;
    for (int mask = 0; mask < (1 << n); mask++){
        if (!(mask >> src & 1) || (mask >> sink & 1)) continue;
        long long cut = 0;
        for (auto& a : arcs) if ((mask >> a.u & 1) && !(mask >> a.v & 1)) cut += a.cap;
        best = min(best, cut);
    }
    return best;
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
        long long flow = g.maxflow();
        assert(flow == min_cut(n, src, sink, arcs));
        check_flow(g, flow);
        assert(dense.maxflow() == flow);
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
    return 0;
}
