#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAX 100010


using namespace std;

typedef pair<int, int> Pair;

struct PairHash{
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1, T2> const &pair) const{
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);
        return h1 ^ h2;
    }
};

struct Bridge{
    int u, v; /// bridge edge from node u to v
    int cnt_u, cnt_v; /// number of nodes in the connected component of u and v if bridge edge is disconnected

    Bridge(){}
    Bridge(int u, int v, int cnt_u, int cnt_v) : u(u), v(v), cnt_u(cnt_u), cnt_v(cnt_v) {}
};


/***
 *
 * finds bridge edges in simple graphs in O(N + M)
 * a bridge is such an edge which when removed makes the graph disconnected
 * or more precisely, increases the number of connected components in the graph
 * nodes are numbered from 0 to N-1
 *
***/

struct Graph{
    bool visited[MAX];
    vector <int> adj[MAX];
    int n, dt, discover[MAX], low[MAX], cmp[MAX], num[MAX];

    Graph() {}
    Graph(int nodes){
        n = nodes;
        for (int i = 0; i < MAX; i++) adj[i].clear();
    }

    void dfs(int u, int p, vector <Bridge> &bridges){
        visited[u] = true;
        discover[u] = low[u] = ++dt;

        for (auto v: adj[u]){
            if (!visited[v]){
                dfs(v, u, bridges);
                low[u] = min(low[u], low[v]);

                if (low[v] > discover[u]){
                    int cnt = dt - discover[v] + 1;
                    bridges.push_back(Bridge(u, v, cmp[u] - cnt, cnt));
                }
            }
            else if (v != p) low[u] = min(low[u], discover[v]);
        }
    }

    void dfs(int u){
        low[dt++] = u;
        visited[u] = true;

        for (auto v: adj[u]){
            if (!visited[v]) dfs(v);
        }
    }

    /// adds undirected edge from u to v
    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector <Bridge> get_bridges(){
        vector <Bridge> bridges;
        memset(visited, 0, sizeof(visited));

        for (int i = 0; i < n; i++){
            if (!visited[i]){
                dt = 0;
                dfs(i);
                for (int j = 0; j < dt; j++) cmp[low[j]] = dt;
            }
        }

        memset(visited, 0, sizeof(visited));
        for (int i = 0; i < n; i++){
            if (!visited[i]){
                dt = 0;
                dfs(i, -1, bridges);
            }
        }

        return bridges;
    }

    void dfs_bridge_tree(int u, int id, tr1::unordered_set <Pair, PairHash>& bridge_set){
        num[u] = id;
        for (auto v: adj[u]) {
            if (num[v] == -1 && !bridge_set.count(Pair(min(u, v), max(u, v)))){
                dfs_bridge_tree(v, id, bridge_set);
            }
        }
    }

    /***
     * removing all the bridges will divide the graph into different bridge components
     * let the number of components in this graph be k
     * label all these components arbitrarily from 0 to k - 1
     * for all nodes in the original graph, let num[x] be the component number for any node x
     * the bridge tree is a new graph formed with the bridge components
     * it is simply the graph where all the edges are (num[u], num[v]) for all (u, v) bridges
     * the number of nodes can be derived implicitly by taking the max node number of all the edges + 1
     *
     * note that although we are calling it a tree, this will actually be a forest of trees
     *
    ***/

    vector <Pair> get_bridge_tree (){
        auto bridges = get_bridges();
        tr1::unordered_set <Pair, PairHash> bridge_set;

        for (auto bridge: bridges){
            int u = bridge.u, v = bridge.v;
            bridge_set.insert(Pair(min(u, v), max(u, v)));
        }

        int id = 0;
        memset(num, -1, sizeof(num));
        for (int u = 0; u < n; u++){
            if (num[u] == -1){
                dfs_bridge_tree(u, id, bridge_set);
                id++;
            }
        }

        vector <Pair> bridge_tree;
        for (auto bridge: bridges){
            int u = bridge.u, v = bridge.v;
            bridge_tree.push_back(Pair(num[u], num[v]));
        }

        return bridge_tree;
    }
};


int main(){
    auto graph = Graph(10);
    graph.add_edge(0, 1);
    graph.add_edge(1, 2);
    graph.add_edge(2, 0);
    graph.add_edge(5, 6);
    graph.add_edge(6, 7);
    graph.add_edge(7, 5);
    graph.add_edge(2, 3);
    graph.add_edge(3, 4);
    graph.add_edge(2, 4);
    graph.add_edge(4, 5);
    graph.add_edge(8, 9);

    auto bridges = graph.get_bridges();
    printf("%d bridge edges found:\n", (int)bridges.size());   /// 2 bridge edges found:
    for (auto bridge: bridges) {
        /// Edge 4-5, disconnecting will decompose the graph into biconnected components of sizes 5 and 3
        /// Edge 8-9, disconnecting will decompose the graph into biconnected components of sizes 1 and 1
        printf("Edge %d-%d, disconnecting will decompose the graph into biconnected components of sizes %d and %d\n", bridge.u, bridge.v, bridge.cnt_u, bridge.cnt_v);
    }
    puts("");

    auto bridge_tree = graph.get_bridge_tree();
    printf("%d edges in the bridge tree:\n", (int)bridge_tree.size());  /// 2 edges in the bridge tree:
    for (auto edge: bridge_tree){
        /// Edge 0-1
        /// Edge 2-3
        printf("Edge %d-%d\n", edge.first, edge.second);
    }

    return 0;
}
