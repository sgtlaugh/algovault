/***
 *
 * Fast maximum flow with Dinic's algorithm
 * 0 based indexing for nodes, so nodes are numbered from 0 to n-1
 *
 * Use FlowGraph for standard flow with edge capacity
 * Use FlowGraphWithNodeCap when nodes can have capacity as well
 * Use DenseFlowGraph for dense graphs: an n x n capacity matrix, same API, O(n^2) memory (32 MB at n = 2000)
 *     measured 2x to 7x faster than FlowGraph with half or more of all n^2 edges present, n from 600 to 2000
 *     same construction and maxflow() calls, but no per-edge flow values since only residual capacities are kept
 *
 * For more speed, get rid of the struct and wrap it up in a namespace or make it global (25% speed gain locally)
 * If you need to initialize the struct many times, make the arrays vectors or get rid of the struct as above
 *
***/

#include <bits/stdc++.h>

#define MAXN 50010

using namespace std;

struct Edge{
    int u, v;
    long long cap, flow;

    Edge(){}
    Edge(int u, int v, long long cap, long long flow) : u(u), v(v), cap(cap), flow(flow) {}
};

struct FlowGraph{
    vector <int> adj[MAXN];
    vector <struct Edge> E;
    int n, src, sink, Q[MAXN], ptr[MAXN], dis[MAXN];

    FlowGraph(){}
    FlowGraph(int n, int src, int sink): n(n), src(src), sink(sink) {}

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
        int u, f = 0, l = 0;
        memset(dis, -1, sizeof(dis[0]) * n);

        dis[src] = 0, Q[l++] = src;
        while (f < l && dis[sink] == -1){
            u = Q[f++];
            for (auto id: adj[u]){
                if (dis[E[id].v] == -1 && E[id].flow < E[id].cap){
                    Q[l++] = E[id].v;
                    dis[E[id].v] = dis[u] + 1;
                }
            }
        }
        return dis[sink] != -1;
    }

    long long dfs(int u, long long flow){
        if (u == sink || !flow) return flow;

        int len = adj[u].size();
        while (ptr[u] < len){
            int id = adj[u][ptr[u]];
            if (dis[E[id].v] == dis[u] + 1){
                long long f = dfs(E[id].v, min(flow, E[id].cap - E[id].flow));
                if (f){
                    E[id].flow += f, E[id ^ 1].flow -= f;
                    return f;
                }
            }
            ptr[u]++;
        }

        return 0;
    }

    long long maxflow(){
        assert(src != sink);
        long long flow = 0;

        while (bfs()){
            memset(ptr, 0, n * sizeof(ptr[0]));
            while (long long f = dfs(src, LLONG_MAX)){
                flow += f;
            }
        }

        return flow;
    }
};

struct FlowGraphWithNodeCap{
    FlowGraph flowgraph;

    FlowGraphWithNodeCap(int n, int src, int sink, vector <long long> node_capacity){
        flowgraph = FlowGraph(2 * n, 2 * src, 2 * sink + 1);

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

    DenseFlowGraph(int n, int src, int sink) : n(n), src(src), sink(sink), cap(n, vector<long long>(n, 0)), dis(n), ptr(n) {}

    void add_directed_edge(int u, int v, long long c){
        cap[u][v] += c;
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
                cap[u][v] -= f, cap[v][u] += f;
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

    assert(flow_graph.maxflow() == 5);

    auto flow_graph_node_cap = FlowGraphWithNodeCap(n, 0, n - 1, {5, 4, 3, 2});

    flow_graph_node_cap.add_edge(0, 1, 3);
    flow_graph_node_cap.add_edge(1, 2, 4);
    flow_graph_node_cap.add_edge(2, 0, 2);
    flow_graph_node_cap.add_edge(1, 1, 5);
    flow_graph_node_cap.add_edge(2, 3, 3);
    flow_graph_node_cap.add_edge(3, 2, 3);

    assert(flow_graph_node_cap.maxflow() == 2);

    return 0;
}
