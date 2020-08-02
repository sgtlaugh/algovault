/***
 *
 * Fast maximum flow with Dinic's algorithm
 * 0 based indexing for nodes, so nodes are numbered from 0 to n-1
 *
 * Use FlowGraph for standard flow with edge capacity
 * Use FlowGraphWithNodeCap when nodes can have capacity as well
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

    void add_edge(int u, int v, int cap){
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

int main(){
    const int n = 4;
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
