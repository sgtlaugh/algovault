/***
 *
 * Fast maximum flow with Dinic's algorithm
 * 0 based indexing for nodes, so nodes are numbered from 0 to n-1
 *
 * Complexity: O(V^2 E), O(E sqrt(V)) on unit-capacity bipartite matching
 *     Library Checker bipartitematching (1e5 + 1e5 nodes, 2e5 edges) slowest case 1.9 s at -O2 against a 5 s limit, on par with hopcroft_karp.cpp
 *     FlowGraph::maxflow() works on flat per-node copies of the edges, about 20 extra bytes per entry of E, and writes the flows back into E
 *     so maxflow() can be called again after adding edges, continuing from the current flow
 *     maxflow() is the only entry point, bfs() and blocking_flow() read those copies and are UB before maxflow() builds them
 *
 * Use FlowGraph for standard flow with edge capacity
 * Use FlowGraphWithNodeCap when nodes can have capacity as well
 * Use DenseFlowGraph for dense graphs: an n x n capacity matrix, same API, O(n^2) memory (32 MB at n = 2000)
 *     measured 2x to 7x faster than FlowGraph with half or more of all n^2 edges present, n from 600 to 2000
 *     same construction and maxflow() calls, but no per-edge flow values since only residual capacities are kept
 *
 * Min cut:
 *     after maxflow(), in_source_side(u) is true for the nodes reachable from src in the residual graph
 *     that is the smallest source side over all min cuts (FlowGraph and DenseFlowGraph)
 *     FlowGraph only: cut_edges() gives the ids into E of the edges leaving that side, their capacities sum to the max flow
 *     FlowGraphWithNodeCap: use .flowgraph.in_source_side() and .flowgraph.cut_edges()
 *         node x is 2x (in) and 2x + 1 (out) there, a cut edge 2x -> 2x + 1 is a cut vertex
 *     before maxflow() the source side is empty and cut_edges() returns nothing
 *
 * Max closure / project selection: pick a set of items with weights w[i], where picking i forces picking j
 *     src -> i with capacity w[i] if w[i] > 0, i -> sink with capacity -w[i] if w[i] < 0, i -> j with capacity LLONG_MAX
 *     best total weight = sum of positive w[i] - maxflow(), the picked items are those with in_source_side(i)
 *     the sum of positive w[i] must be below LLONG_MAX, the LLONG_MAX edges themselves never overflow
 *     DenseFlowGraph caps each merged u <-> v cell at LLONG_MAX, exact while the max flow is below LLONG_MAX
 *
 * Max density subgraph, max |E(S)| / |S| over nonempty S: Dinkelbach on the closure above, n + m + 2 nodes
 *     start with the density a / b = m / n of the whole graph
 *     closure items: each edge with weight b, forcing both endpoints, each vertex with weight -a
 *     if the best closure weight is positive, set a / b to |E(S)| / |S| for its picked vertices S and repeat, else a / b is the answer
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Edge{
    int u, v;
    long long cap, flow;

    Edge() {}
    Edge(int u, int v, long long cap, long long flow) : u(u), v(v), cap(cap), flow(flow) {}
};

struct FlowGraph{
    int n, src, sink;
    vector<vector<int>> adj;
    vector<struct Edge> E;
    vector<int> Q, ptr, dis, path, start, to, rev, edge_id;
    vector<long long> res;

    FlowGraph(int n, int src, int sink): n(n), src(src), sink(sink), adj(n), Q(n), ptr(n), dis(n, -1) {}

    void add_directed_edge(int u, int v, long long cap){
        adj[u].push_back(E.size());
        E.push_back(Edge(u, v, cap, 0));

        adj[v].push_back(E.size());
        E.push_back(Edge(v, u, 0, 0));
    }

    void add_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    bool bfs(){
        int f = 0, l = 0;
        fill(dis.begin(), dis.end(), -1);

        dis[src] = 0, Q[l++] = src;
        while (f < l && dis[sink] == -1){
            int u = Q[f++];
            for (int i = start[u]; i < start[u + 1]; i++){
                if (res[i] > 0 && dis[to[i]] == -1) dis[to[i]] = dis[u] + 1, Q[l++] = to[i];
            }
        }
        return dis[sink] != -1;
    }

    /// Walks from sink back to src so it only enters nodes the bfs reached from src, a walk from src spent most of its time in branches
    /// that never reach sink, 2.5x to 3x slower at -O2 on Library Checker's cycle and unique matching cases. Iterative so that a src to sink distance of n cannot overflow the stack
    /// path[j] is the arc scanned at the j-th node from sink, flow goes along its reverse rev[path[j]]
    long long blocking_flow(){
        long long flow = 0;
        int u = sink;
        path.clear();

        while (true){
            if (u == src){
                long long f = LLONG_MAX;
                for (int i : path) f = min(f, res[rev[i]]);

                int k = path.size();
                for (int j = k - 1; j >= 0; j--){
                    res[rev[path[j]]] -= f, res[path[j]] += f;
                    if (res[rev[path[j]]] == 0) k = j;
                }
                flow += f, u = to[rev[path[k]]];
                path.resize(k);
                continue;
            }

            while (ptr[u] < start[u + 1] && (res[rev[ptr[u]]] == 0 || dis[to[ptr[u]]] != dis[u] - 1)) ptr[u]++;

            if (ptr[u] < start[u + 1]){
                path.push_back(ptr[u]);
                u = to[ptr[u]];
            }
            else{
                if (u == sink) return flow;
                u = to[rev[path.back()]];
                path.pop_back();
                ptr[u]++;
            }
        }
    }

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 with residual capacity res, reading E through adj ids missed the cache on every arc
    void build_residual(){
        int m = E.size();
        vector<int> pos(m);
        start.assign(n + 1, m), to.resize(m), rev.resize(m), edge_id.resize(m), res.resize(m);

        for (int u = 0, i = 0; u < n; u++){
            start[u] = i;
            for (int id : adj[u]) pos[id] = i, edge_id[i] = id, to[i] = E[id].v, res[i] = E[id].cap - E[id].flow, i++;
        }
        for (int i = 0; i < m; i++) rev[i] = pos[edge_id[i] ^ 1];
    }

    long long maxflow(){
        assert(src != sink);
        long long flow = 0;
        build_residual();

        while (bfs()){
            copy(start.begin(), start.end() - 1, ptr.begin());
            flow += blocking_flow();
        }

        for (int i = 0; i < (int)res.size(); i++) E[edge_id[i]].flow = E[edge_id[i]].cap - res[i];
        return flow;
    }

    /// Reads dis from the final bfs of maxflow(), which failed to reach sink and so visited every residual-reachable node
    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    vector<int> cut_edges() const{
        vector<int> ids;
        for (int id = 0; id < (int)E.size(); id++){
            if (E[id].cap > 0 && in_source_side(E[id].u) && !in_source_side(E[id].v)) ids.push_back(id);
        }
        return ids;
    }
};

struct FlowGraphWithNodeCap{
    FlowGraph flowgraph;

    FlowGraphWithNodeCap(int n, int src, int sink, vector<long long> node_capacity) : flowgraph(2 * n, 2 * src, 2 * sink + 1){
        for (int i = 0; i < n; i++){
            flowgraph.add_directed_edge(2 * i, 2 * i + 1, node_capacity[i]);
        }
    }

    void add_directed_edge(int u, int v, long long cap){
        flowgraph.add_directed_edge(2 * u + 1, 2 * v, cap);
    }

    void add_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(){
        return flowgraph.maxflow();
    }
};

/// Dinic on an adjacency matrix: no edge lists, the residual capacity of u -> v is cap[u][v]
struct DenseFlowGraph{
    int n, src, sink;
    vector<vector<long long>> cap;
    vector<int> dis, ptr;

    DenseFlowGraph(int n, int src, int sink) : n(n), src(src), sink(sink), cap(n, vector<long long>(n, 0)), dis(n, -1), ptr(n) {}

    void add_directed_edge(int u, int v, long long c){
        saturating_add(cap[u][v], c);
    }

    void add_edge(int u, int v, long long c){
        add_directed_edge(u, v, c);
        add_directed_edge(v, u, c);
    }

    bool bfs(){
        fill(dis.begin(), dis.end(), -1);
        vector<int> queue = {src};
        dis[src] = 0;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (int v = 0; v < n; v++){
                if (dis[v] == -1 && cap[u][v] > 0) dis[v] = dis[u] + 1, queue.push_back(v);
            }
        }
        return dis[sink] != -1;
    }

    long long dfs(int u, long long flow){
        if (u == sink || !flow) return flow;
        for (int& v = ptr[u]; v < n; v++){
            if (dis[v] != dis[u] + 1 || cap[u][v] <= 0) continue;
            long long f = dfs(v, min(flow, cap[u][v]));
            if (f){
                cap[u][v] -= f, saturating_add(cap[v][u], f);
                return f;
            }
        }
        return 0;
    }

    long long maxflow(){
        assert(src != sink);
        long long flow = 0;

        while (bfs()){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, LLONG_MAX)) flow += f;
        }
        return flow;
    }

    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    /// One cell holds u -> v plus the flow on v -> u, so two LLONG_MAX edges would overflow it, LLONG_MAX still acts as infinite
    static void saturating_add(long long& cell, long long c){
        cell = cell > LLONG_MAX - c ? LLONG_MAX : cell + c;
    }
};

int main(){
    const int n = 4;
    auto dense = DenseFlowGraph(n, 0, n - 1);
    dense.add_edge(0, 1, 3), dense.add_edge(1, 2, 4), dense.add_edge(2, 0, 2), dense.add_edge(1, 1, 5), dense.add_edge(2, 3, 3), dense.add_edge(3, 2, 3);
    assert(dense.maxflow() == 5);

    auto flow_graph = FlowGraph(n, 0, n - 1);

    flow_graph.add_edge(0, 1, 3);
    flow_graph.add_edge(1, 2, 4);
    flow_graph.add_edge(2, 0, 2);
    flow_graph.add_edge(1, 1, 5);
    flow_graph.add_edge(2, 3, 3);
    flow_graph.add_edge(3, 2, 3);

    assert(flow_graph.cut_edges().empty());
    for (int u = 0; u < n; u++) assert(!flow_graph.in_source_side(u) && !DenseFlowGraph(n, 0, n - 1).in_source_side(u));

    assert(flow_graph.maxflow() == 5);
    for (int u = 0; u < n; u++) assert(flow_graph.in_source_side(u) == (u == 0) && dense.in_source_side(u) == (u == 0));

    set<tuple<int, int, long long>> cut;
    for (int id : flow_graph.cut_edges()) cut.insert({flow_graph.E[id].u, flow_graph.E[id].v, flow_graph.E[id].cap});
    assert((cut == set<tuple<int, int, long long>>{{0, 1, 3}, {0, 2, 2}}));

    flow_graph.add_directed_edge(0, 3, 4);
    assert(flow_graph.maxflow() == 4 && flow_graph.E.back().flow == -4);

    auto flow_graph_node_cap = FlowGraphWithNodeCap(n, 0, n - 1, {5, 4, 3, 2});

    flow_graph_node_cap.add_edge(0, 1, 3);
    flow_graph_node_cap.add_edge(1, 2, 4);
    flow_graph_node_cap.add_edge(2, 0, 2);
    flow_graph_node_cap.add_edge(1, 1, 5);
    flow_graph_node_cap.add_edge(2, 3, 3);
    flow_graph_node_cap.add_edge(3, 2, 3);

    assert(flow_graph_node_cap.maxflow() == 2);
    auto node_cut = flow_graph_node_cap.flowgraph.cut_edges();
    assert(node_cut.size() == 1 && flow_graph_node_cap.flowgraph.E[node_cut[0]].u == 6 && flow_graph_node_cap.flowgraph.E[node_cut[0]].v == 7);

    vector<long long> weight = {10, 5, -4, -3, -6};
    auto closure = FlowGraph(7, 5, 6);
    for (int i = 0; i < 5; i++){
        if (weight[i] > 0) closure.add_directed_edge(5, i, weight[i]);
        else closure.add_directed_edge(i, 6, -weight[i]);
    }
    closure.add_directed_edge(0, 2, LLONG_MAX), closure.add_directed_edge(0, 3, LLONG_MAX);
    closure.add_directed_edge(1, 3, LLONG_MAX), closure.add_directed_edge(1, 4, LLONG_MAX);

    assert(15 - closure.maxflow() == 3);
    for (int i = 0; i < 5; i++) assert(closure.in_source_side(i) == (i == 0 || i == 2 || i == 3));

    return 0;
}
