/***
 *
 * Mincost Maxflow with SPFA
 *
 * Initialize with the maximum number of nodes in the graph first
 * Works both for 0 and 1 based indexing
 *
 * Complexity is exponential in theory, however in practice it is very fast
 * Further reading - https://codeforces.com/blog/entry/70740
 *
 * Note, flow[2*i] contains the amount of flow in i-th edge
 * So if we ever need to print it we can retrieve it easily
 *
***/

#include <bits/stdc++.h>

using namespace std;

const int MAXV = 100010;
const int MAXE = 1000010;

template <typename T>
struct FlowGraph{
    const T INF = numeric_limits<T>::max();

    queue<int> Q;
    bitset<MAXV> visited;
    int n, id, head[MAXV], nxt[MAXE], to[MAXE];
    T cap[MAXE], cost[MAXE], flow[MAXE], dis[MAXV];

    FlowGraph(int n) : n(n), id(1){
        for (int i = 0; i <= n; i++) head[i] = nxt[i] = 0;
    }

    /// Add an edge from u to v with edge capacity c and edge weight w
    void add_edge(int u, int v, T c, T w, bool directed=true){
        nxt[++id] = head[u], to[id] = v, head[u] = id, cap[id] = c, cost[id] = w, flow[id] = 0;
        nxt[++id] = head[v], to[id] = u, head[v] = id, cap[id] = 0, cost[id] = -w, flow[id] = 0;
        if (!directed) add_edge(v, u, v, w, true);
    }

    bool spfa(int src, int sink){
        int i, u, v;
        for (i = 0; i <= n; i++) dis[i] = INF, visited[i] = 0;

        dis[src] = 0;
        if (!visited[src]) Q.push(src), visited[src] = 1;

        while (!Q.empty()){
            u = Q.front();
            Q.pop(), visited[u] = 0;

            for (i = head[u]; i; i = nxt[i]){
                v = to[i];
                auto new_dis = dis[u] + cost[i];

                if (flow[i] < cap[i] && new_dis < dis[v]){
                    dis[v] = new_dis;
                    if (!visited[v]){
                        Q.push(v), visited[v] = 1;
                    }
                }
            }
        }

        return dis[sink] < INF;
    }

    T dfs(int u, int sink, T lim){
        visited[u] = 1;
        if (u == sink) return lim;

        T aug = 0;
        for (int i = head[u]; i && aug < lim; i = nxt[i]){
            int v = to[i];
            if (!visited[v] && flow[i] < cap[i] && dis[v] == dis[u] + cost[i]){
                auto x = dfs(v, sink, min(cap[i] - flow[i], lim - aug));
                flow[i] += x, flow[i ^ 1] -= x;
                aug += x;
            }
        }

        return aug;
    }

    pair<T, T> mincost_maxflow(int src, int sink){
        T mincost = 0, maxflow = 0;

        while (spfa(src, sink)){
            T aug = INF;
            while (aug){
                for (int i = 0; i <= n; i++) visited[i] = 0;
                aug = dfs(src, sink, INF);
                maxflow += aug, mincost += aug * dis[sink];
            }
        }

        return {mincost, maxflow};
    }
};

int main(){
    const int n = 6;
    auto flow_graph = FlowGraph<int>(n);

    /// Graph taken from https://cp-algorithms.com/graph/edmonds_karp.html
    flow_graph.add_edge(1, 2, 7, 1);
    flow_graph.add_edge(2, 3, 5, 2);
    flow_graph.add_edge(3, 6, 8, 3);
    flow_graph.add_edge(1, 4, 4, 4);
    flow_graph.add_edge(4, 2, 3, 5);
    flow_graph.add_edge(4, 5, 2, 4);
    flow_graph.add_edge(2, 5, 3, 3);
    flow_graph.add_edge(5, 3, 3, 2);
    flow_graph.add_edge(5, 6, 5, 1);

    int mincost, maxflow;
    tie(mincost, maxflow) = flow_graph.mincost_maxflow(1, n);

    assert(mincost == 71);
    assert(maxflow == 10);

    return 0;
}
