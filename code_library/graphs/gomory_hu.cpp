/***
 *
 * Gomory-Hu Tree (Gusfield's algorithm)
 * All-pairs minimum cut of an undirected graph from n - 1 max flows, 0-based nodes
 *
 * Complexity: n - 1 Dinic max flows, each O(n^2 m) worst case, O(n + m) memory
 *             O(n) per min_cut query
 *
 * GomoryHu g(n); g.add_edge(u, v, cap); ... g.build();
 * g.min_cut(u, v): minimum cut between u != v, the smallest weight on the tree path between them
 * After build, the tree has edges (v, tree_parent[v]) of weight tree_weight[v] for v = 1..n-1, rooted at node 0
 *     removing edge (v, tree_parent[v]) leaves the subtree of v on one side, a minimum cut between v and
 *     tree_parent[v] in the graph of value tree_weight[v], so subtrees give min cut partitions, not only values
 * Parallel edges and self loops are allowed, a disconnected pair has min cut 0
 * Total capacity must fit in long long, Dinic's dfs recurses up to n deep
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct GomoryHu{
    struct Edge{
        int v;
        long long cap, flow;
    };

    int n;
    vector<vector<int>> adj;
    vector<Edge> E;
    vector<int> Q, ptr, dis, tree_parent, depth;
    vector<long long> tree_weight;

    GomoryHu(int n) : n(n), adj(n), Q(n), ptr(n), dis(n) {}

    /// One edge pair with capacity cap both ways: pushing f along one side frees f on the other
    void add_edge(int u, int v, long long cap){
        adj[u].push_back(E.size());
        E.push_back({v, cap, 0});

        adj[v].push_back(E.size());
        E.push_back({u, cap, 0});
    }

    void build(){
        tree_parent.assign(n, 0);
        tree_weight.assign(n, 0);

        for (int s = 1; s < n; s++){
            int t = tree_parent[s];
            for (auto& e : E) e.flow = 0;
            tree_weight[s] = maxflow(s, t);

            for (int v = 0; v < n; v++){
                if (v != s && tree_parent[v] == t && dis[v] != -1) tree_parent[v] = s;
            }

            /// Gusfield's swap: without it the tree is only flow-equivalent and subtrees are not min cuts
            if (dis[tree_parent[t]] != -1){
                tree_parent[s] = tree_parent[t], tree_parent[t] = s;
                swap(tree_weight[s], tree_weight[t]);
            }
        }

        vector<vector<int>> children(n);
        for (int v = 1; v < n; v++) children[tree_parent[v]].push_back(v);
        depth.assign(n, 0);
        vector<int> order = {0};
        for (int i = 0; i < (int)order.size(); i++){
            for (int c : children[order[i]]) depth[c] = depth[order[i]] + 1, order.push_back(c);
        }
    }

    long long min_cut(int u, int v) const{
        assert(u != v && (int)depth.size() == n);
        long long res = LLONG_MAX;
        while (u != v){
            if (depth[u] < depth[v]) swap(u, v);
            res = min(res, tree_weight[u]);
            u = tree_parent[u];
        }
        return res;
    }

    bool bfs(int src, int sink){
        int f = 0, l = 0;
        fill(dis.begin(), dis.end(), -1);

        dis[src] = 0, Q[l++] = src;
        while (f < l && dis[sink] == -1){
            int u = Q[f++];
            for (auto id : adj[u]){
                if (dis[E[id].v] == -1 && E[id].flow < E[id].cap){
                    Q[l++] = E[id].v;
                    dis[E[id].v] = dis[u] + 1;
                }
            }
        }
        return dis[sink] != -1;
    }

    long long dfs(int u, int sink, long long flow){
        if (u == sink || !flow) return flow;

        int len = adj[u].size();
        while (ptr[u] < len){
            int id = adj[u][ptr[u]];
            if (dis[E[id].v] == dis[u] + 1){
                long long f = dfs(E[id].v, sink, min(flow, E[id].cap - E[id].flow));
                if (f){
                    E[id].flow += f, E[id ^ 1].flow -= f;
                    return f;
                }
            }
            ptr[u]++;
        }

        return 0;
    }

    /// On return dis[v] != -1 exactly for the source side of a minimum cut, as the last bfs failed
    long long maxflow(int src, int sink){
        long long flow = 0;

        while (bfs(src, sink)){
            fill(ptr.begin(), ptr.end(), 0);
            while (long long f = dfs(src, sink, LLONG_MAX)) flow += f;
        }

        return flow;
    }
};

int main(){
    /***
     *           3
     *       0 ----- 1
     *        \     /
     *       1 \   / 2
     *          \ /
     *           2 ----- 3
     *               4
    ***/
    GomoryHu g(4);
    g.add_edge(0, 1, 3);
    g.add_edge(1, 2, 2);
    g.add_edge(0, 2, 1);
    g.add_edge(2, 3, 4);
    g.build();
    assert(g.min_cut(0, 1) == 4);  /// cut off 0: 3 + 1
    assert(g.min_cut(1, 3) == 3);  /// {0, 1} against {2, 3}: 1 + 2
    assert(g.min_cut(2, 3) == 4);

    /// Each tree edge (v, tree_parent[v]) carries the min cut between its endpoints
    for (int v = 1; v < 4; v++) assert(g.tree_weight[v] == g.min_cut(v, g.tree_parent[v]));

    GomoryHu split(4);
    split.add_edge(0, 1, 5);
    split.add_edge(2, 3, 7);
    split.build();
    assert(split.min_cut(0, 2) == 0);  /// different components
    return 0;
}
