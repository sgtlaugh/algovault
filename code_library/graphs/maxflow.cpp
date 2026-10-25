/***
 *
 * Maximum flow with Dinic's algorithm
 * Max flow and min cut on nodes 0 to n - 1, with optional node capacities or a dense capacity matrix
 *
 * Complexity: O(n^2 m) for m edges, O(m sqrt(n)) on unit-capacity bipartite matching
 *
 * FlowGraph g(n, src, sink): empty graph with edge capacities, src != sink
 * g.add_directed_edge(u, v, cap): adds u -> v and returns its id into g.E
 * g.add_undirected_edge(u, v, cap): adds u -> v and v -> u, both with capacity cap
 * g.maxflow(limit = LLONG_MAX): pushes as much flow as possible up to limit and returns the flow this call added
 * g.cut_edges(): ids into E of the edges leaving the source side, their capacities sum to the max flow
 * g.in_source_side(u): after maxflow(), true if src reaches u in the residual graph, the smallest min cut side
 * g.E[id].flow: flow on edge id after maxflow()
 *
 * FlowGraphWithNodeCap g(n, src, sink, node_capacity): the edge and maxflow() calls of FlowGraph, on g.flowgraph
 * DenseFlowGraph g(n, src, sink): an n x n capacity matrix, FlowGraph's calls without cut_edges() or per-edge flows
 *
 * Notes:
 *   - On bipartite matching hopcroft_karp.cpp is up to 2x faster on random graphs, on par on Library Checker's slowest
 *   - maxflow() keeps flat per-node copies of the edges, about 20 extra bytes per edge
 *   - maxflow() is the only entry point, bfs() and blocking_flow() are UB before it builds the flat copies
 *   - Before maxflow() the source side is empty and cut_edges() returns nothing
 *   - A call stopped by its limit below the max flow leaves in_source_side(sink) true and cut_edges() no min cut
 *   - Each maxflow() call continues from the flow already pushed, so an unlimited call after a limited one completes it
 *   - Node capacities apply to src and sink too, give them a huge capacity to leave them uncapped
 *   - FlowGraphWithNodeCap edge ids index g.flowgraph.E, where node x is 2x (in) and 2x + 1 (out)
 *   - Its min cut comes from g.flowgraph, a cut edge 2x -> 2x + 1 is a cut vertex x
 *   - Use DenseFlowGraph with half or more of all n^2 edges: 2x to 7x faster for n = 600 to 2000, 32 MB at n = 2000
 *   - Max closure, items with weights w[i] where picking i forces picking j:
 *       - src -> i with capacity w[i] if w[i] > 0, i -> sink with capacity -w[i] if w[i] < 0, i -> j with LLONG_MAX
 *       - The best total weight is the sum of positive w[i] - maxflow(), the picked items are those in_source_side(i)
 *       - The sum of positive w[i] must be below LLONG_MAX, the LLONG_MAX edges never overflow, in DenseFlowGraph too
 *   - Max density subgraph, max |E(S)| / |S| over nonempty S: Dinkelbach on the closure above, n + m + 2 nodes
 *       - Start with the density a / b = m / n of the whole graph
 *       - Closure items: each edge with weight b, forcing both endpoints, each vertex with weight -a
 *       - If the best closure is positive, set a / b to |E(S)| / |S| of its vertices S and repeat, else stop
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FlowGraph{
    struct Edge{
        int u, v;
        long long cap, flow;
    };

    int n, src, sink;
    vector<vector<int>> adj;
    vector<Edge> E;
    vector<int> Q, ptr, dis, path, start, to, rev, edge_id;
    vector<long long> res;

    FlowGraph(int n, int src, int sink): n(n), src(src), sink(sink), adj(n), Q(n), ptr(n), dis(n, -1) {}

    int add_directed_edge(int u, int v, long long cap){
        adj[u].push_back(E.size()), E.push_back({u, v, cap, 0});
        adj[v].push_back(E.size()), E.push_back({v, u, 0, 0});
        return E.size() - 2;
    }

    void add_undirected_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(long long limit = LLONG_MAX){
        assert(src != sink);
        long long flow = 0;
        build_residual();

        /// bfs() runs before the limit check, so in_source_side(sink) tells whether the flow is maximum
        while (bfs() && flow < limit){
            copy(start.begin(), start.end() - 1, ptr.begin());
            flow += blocking_flow(limit - flow);
        }

        for (int i = 0; i < (int)res.size(); i++) E[edge_id[i]].flow = E[edge_id[i]].cap - res[i];
        return flow;
    }

    vector<int> cut_edges() const{
        vector<int> ids;
        for (int id = 0; id < (int)E.size(); id++){
            if (E[id].cap > 0 && in_source_side(E[id].u) && !in_source_side(E[id].v)) ids.push_back(id);
        }
        return ids;
    }

    /// Reads dis from the final bfs of maxflow(), which failed to reach sink and so visited every node src reaches
    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    /// Arcs of u sit at start[u] .. start[u + 1] - 1 with residual res: reading E through adj missed the cache
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

    /// Iterative walk from sink enters only nodes src reaches (2.5x-3x faster), path[j] is the arc taken j nodes back
    long long blocking_flow(long long limit){
        long long flow = 0;
        int u = sink;

        path.clear();

        while (true){
            if (u == src){
                long long f = limit - flow;
                for (int i : path) f = min(f, res[rev[i]]);

                int k = path.size();
                for (int j = k - 1; j >= 0; j--){
                    res[rev[path[j]]] -= f, res[path[j]] += f;
                    if (res[rev[path[j]]] == 0) k = j;
                }
                flow += f;
                if (flow == limit) return flow;

                u = to[rev[path[k]]];
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
};

struct FlowGraphWithNodeCap{
    FlowGraph flowgraph;

    FlowGraphWithNodeCap(int n, int src, int sink, vector<long long> node_capacity) : flowgraph(2 * n, 2 * src, 2 * sink + 1){
        assert(src != sink);
        for (int i = 0; i < n; i++){
            flowgraph.add_directed_edge(2 * i, 2 * i + 1, node_capacity[i]);
        }
    }

    int add_directed_edge(int u, int v, long long cap){
        return flowgraph.add_directed_edge(2 * u + 1, 2 * v, cap);
    }

    void add_undirected_edge(int u, int v, long long cap){
        add_directed_edge(u, v, cap);
        add_directed_edge(v, u, cap);
    }

    long long maxflow(long long limit = LLONG_MAX){
        return flowgraph.maxflow(limit);
    }
};

/// cap[u][v] is the residual capacity of u -> v, no flow is stored
struct DenseFlowGraph{
    int n, src, sink;
    vector<vector<long long>> cap;
    vector<int> dis, ptr;

    DenseFlowGraph(int n, int src, int sink) : n(n), src(src), sink(sink), cap(n, vector<long long>(n, 0)), dis(n, -1), ptr(n) {}

    void add_directed_edge(int u, int v, long long c){
        saturating_add(cap[u][v], c);
    }

    void add_undirected_edge(int u, int v, long long c){
        add_directed_edge(u, v, c);
        add_directed_edge(v, u, c);
    }

    long long maxflow(long long limit = LLONG_MAX){
        assert(src != sink);
        long long flow = 0;

        while (bfs() && flow < limit){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, limit - flow)) flow += f;
        }

        return flow;
    }

    bool in_source_side(int u) const{
        return dis[u] != -1;
    }

    /// One cell holds u -> v plus the flow on v -> u, so it saturates at LLONG_MAX, which still means infinite
    static void saturating_add(long long& cell, long long c){
        cell = cell > LLONG_MAX - c ? LLONG_MAX : cell + c;
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
};

int main(){
    /// Max flow 5 from 0 to 3, the two edges out of 0 are the min cut
    FlowGraph g(4, 0, 3);
    int first = g.add_directed_edge(0, 1, 3);
    int second = g.add_directed_edge(0, 2, 2);
    g.add_directed_edge(1, 2, 4);
    g.add_directed_edge(1, 3, 2);
    g.add_directed_edge(2, 3, 4);
    assert(g.maxflow() == 5);
    assert(g.E[first].flow == 3);
    assert((g.cut_edges() == vector<int>{first, second}));
    assert(g.in_source_side(0));
    assert(!g.in_source_side(1));

    int direct = g.add_directed_edge(0, 3, 4);
    assert(g.maxflow() == 4);
    assert(g.E[direct].flow == 4);

    FlowGraphWithNodeCap capped(4, 0, 3, {10, 1, 10, 10});
    capped.add_directed_edge(0, 1, 3);
    int bypass = capped.add_directed_edge(0, 2, 2);
    capped.add_directed_edge(1, 2, 4);
    capped.add_directed_edge(1, 3, 2);
    capped.add_directed_edge(2, 3, 4);
    assert(capped.maxflow() == 3);
    /// Edge 2 is node 1's split edge 2 -> 3, so node 1 is a cut vertex
    assert((capped.flowgraph.cut_edges() == vector<int>{2, bypass}));

    /// The first call stops at its limit of 2, below the max flow 3, so it leaves no cut and the next call adds 1
    FlowGraph limited(3, 0, 2);
    limited.add_directed_edge(0, 1, 5);
    limited.add_directed_edge(1, 2, 3);
    assert(limited.maxflow(2) == 2);
    assert(limited.in_source_side(2));
    assert(limited.maxflow() == 1);
    assert(!limited.in_source_side(2));

    DenseFlowGraph dense(4, 0, 3);
    dense.add_directed_edge(0, 1, 3);
    dense.add_directed_edge(0, 2, 2);
    dense.add_directed_edge(1, 2, 4);
    dense.add_directed_edge(1, 3, 2);
    dense.add_directed_edge(2, 3, 4);
    assert(dense.maxflow() == 5);
    assert(!dense.in_source_side(1));

    /// Max closure: item 0 forces 2 and 3, item 1 forces 3 and 4, so picking {0, 2, 3} is best at 10 - 4 - 3
    vector<long long> weight = {10, 5, -4, -3, -6};
    FlowGraph closure(7, 5, 6);
    for (int i = 0; i < 5; i++){
        if (weight[i] > 0) closure.add_directed_edge(5, i, weight[i]);
        else closure.add_directed_edge(i, 6, -weight[i]);
    }
    closure.add_directed_edge(0, 2, LLONG_MAX);
    closure.add_directed_edge(0, 3, LLONG_MAX);
    closure.add_directed_edge(1, 3, LLONG_MAX);
    closure.add_directed_edge(1, 4, LLONG_MAX);
    assert(15 - closure.maxflow() == 3);
    for (int i = 0; i < 5; i++) assert(closure.in_source_side(i) == (i == 0 || i == 2 || i == 3));

    return 0;
}
