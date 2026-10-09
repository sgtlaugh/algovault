/***
 *
 * Finds bridge edges and builds the bridge tree in undirected graphs, parallel edges and self loops allowed
 * A doubled edge is never a bridge: the search skips only the edge it arrived by, not every edge to the parent
 * add_edge returns the edge index, Bridge::id reports it, so weights or labels can be kept outside
 *
 * A bridge is such an edge which when removed makes the graph disconnected
 * Or more precisely, increases the number of connected components
 *
 * Graph nodes are numbered from 0 to N-1
 *
 * Complexity: O(N + M)
 *
***/

#include <bits/stdtr1c++.h>

#define MAX 100010

using namespace std;

typedef pair<int, int> Pair;

struct Bridge{
    int u, v; /// bridge edge from node u to v
    int cnt_u, cnt_v; /// number of nodes in the connected component of u and v if bridge edge is disconnected
    int id; /// index of the edge, in the order add_edge was called

    Bridge(){}
    Bridge(int u, int v, int cnt_u, int cnt_v, int id = -1) : u(u), v(v), cnt_u(cnt_u), cnt_v(cnt_v), id(id) {}
};

struct Graph{
    bool visited[MAX];
    vector <Pair> adj[MAX]; /// (neighbor, edge index)
    int n, m = 0, dt, discover[MAX], low[MAX], num[MAX];

    Graph() {}
    Graph(int n): n(n) {}

    void dfs(int u, int parent_edge, vector <Bridge> &bridges){
        visited[u] = true;
        discover[u] = low[u] = ++dt;

        for (auto [v, id]: adj[u]){
            if (!visited[v]){
                dfs(v, id, bridges);
                low[u] = min(low[u], low[v]);

                /// emplace_back, push_back(Bridge(...)) grows the GCC -O2 frame enough to overflow 8 MB on a 1e5-node path
                if (low[v] > discover[u]) bridges.emplace_back(u, v, 0, dt - discover[v] + 1, id);
            }
            else if (id != parent_edge) low[u] = min(low[u], discover[v]);
        }
    }

    /// adds undirected edge from u to v and returns its index
    int add_edge(int u, int v){
        adj[u].push_back(Pair(v, m));
        adj[v].push_back(Pair(u, m));
        return m++;
    }

    vector <Bridge> get_bridges(){
        vector <Bridge> bridges;
        memset(visited, 0, sizeof(visited));

        for (int i = 0; i < n; i++){
            if (!visited[i]){
                int first = bridges.size();
                dt = 0;
                dfs(i, -1, bridges);
                /// The component size is known only once its DFS ends
                for (int j = first; j < (int)bridges.size(); j++) bridges[j].cnt_u = dt - bridges[j].cnt_v;
            }
        }

        return bridges;
    }

    void dfs_bridge_tree(int u, int label, const vector <char>& is_bridge){
        num[u] = label;
        for (auto [v, id]: adj[u]) {
            if (num[v] == -1 && !is_bridge[id]){
                dfs_bridge_tree(v, label, is_bridge);
            }
        }
    }

    /***
     * https://www.quora.com/q/threadsiiithyderabad/The-Bridge-Tree-of-a-graph
     *
     * Removing all the bridges will divide the graph into different bridge components
     * Let the number of components in this graph be k
     * Label all these components arbitrarily from 0 to k - 1
     * For all nodes in the original graph, let num[x] be the component number for any node x
     * The bridge tree is a new graph formed with the bridge components
     * It is simply the graph where all the edges are (num[u], num[v]) for all (u, v) bridges
     * The number of nodes can be derived implicitly by taking the max node number of all the edges + 1
     *
     * Note that although we are calling it a tree, this will actually be a forest of trees
     *
    ***/

    vector <Pair> get_bridge_tree (){
        auto bridges = get_bridges();
        vector <char> is_bridge(m, 0);
        for (auto bridge: bridges) is_bridge[bridge.id] = 1;

        int label = 0;
        memset(num, -1, sizeof(num));
        for (int u = 0; u < n; u++){
            if (num[u] == -1){
                dfs_bridge_tree(u, label, is_bridge);
                label++;
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
    /// Each Graph holds about 4 MB of arrays, too much for the stack
    static Graph graph(10);

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
    assert((int)bridges.size() == 2);
    assert(bridges[0].u == 4 && bridges[0].v == 5 && bridges[0].cnt_u == 5 && bridges[0].cnt_v == 3);
    assert(bridges[1].u == 8 && bridges[1].v == 9 && bridges[1].cnt_u == 1 && bridges[1].cnt_v == 1);

    auto bridge_tree = graph.get_bridge_tree();
    assert((int)bridge_tree.size() == 2);
    assert(bridge_tree[0] == Pair(0, 1));
    assert(bridge_tree[1] == Pair(2, 3));
    assert(bridges[0].id == 9 && bridges[1].id == 10);

    static Graph multi(4);
    multi.add_edge(0, 1);
    multi.add_edge(0, 1);
    multi.add_edge(1, 2);
    multi.add_edge(2, 2);
    multi.add_edge(2, 3);
    auto multi_bridges = multi.get_bridges();
    assert((int)multi_bridges.size() == 2);
    assert(multi_bridges[0].id == 4 && multi_bridges[1].id == 2);
    assert(multi_bridges[0].cnt_u == 3 && multi_bridges[0].cnt_v == 1);

    return 0;
}
