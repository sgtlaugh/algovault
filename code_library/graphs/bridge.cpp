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
 * Both searches are iterative with heap-allocated stacks, so a 2e5-node path needs no extra call stack
 *
 * Complexity: O(N + M)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

typedef pair<int, int> Pair;

struct Bridge{
    int u, v; /// bridge edge from node u to v
    int cnt_u, cnt_v; /// number of nodes in the connected component of u and v if bridge edge is disconnected
    int id; /// index of the edge, in the order add_edge was called

    Bridge() {}
    Bridge(int u, int v, int cnt_u, int cnt_v, int id = -1) : u(u), v(v), cnt_u(cnt_u), cnt_v(cnt_v), id(id) {}
};

struct Graph{
    int n, m = 0;
    vector<vector<Pair>> adj; /// (neighbor, edge index)
    vector<int> num; /// bridge component of each node, filled by get_bridge_tree

    Graph(int n): n(n), adj(n), num(n, -1) {}

    /// adds undirected edge from u to v and returns its index
    int add_edge(int u, int v){
        adj[u].push_back(Pair(v, m));
        adj[v].push_back(Pair(u, m));
        return m++;
    }

    /// Bridges are listed in the order the DFS finishes their endpoint v
    vector<Bridge> get_bridges(){
        vector<Bridge> bridges;
        vector<int> discover(n, 0), low(n), parent_edge(n), edge_pos(n, 0), dfs_stack;

        for (int i = 0; i < n; i++){
            if (discover[i]) continue;

            int first = bridges.size(), dt = 0;
            discover[i] = low[i] = ++dt, parent_edge[i] = -1;
            dfs_stack.push_back(i);

            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                if (edge_pos[u] < (int)adj[u].size()){
                    auto [v, id] = adj[u][edge_pos[u]++];
                    if (!discover[v]){
                        discover[v] = low[v] = ++dt, parent_edge[v] = id;
                        dfs_stack.push_back(v);
                    }
                    else if (id != parent_edge[u]) low[u] = min(low[u], discover[v]);
                    continue;
                }

                dfs_stack.pop_back();
                if (dfs_stack.empty()) break;

                int p = dfs_stack.back();
                low[p] = min(low[p], low[u]);
                if (low[u] > discover[p]) bridges.emplace_back(p, u, 0, dt - discover[u] + 1, parent_edge[u]);
            }

            /// The component size is known only once its DFS ends
            for (int j = first; j < (int)bridges.size(); j++) bridges[j].cnt_u = dt - bridges[j].cnt_v;
        }

        return bridges;
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

    vector<Pair> get_bridge_tree (){
        auto bridges = get_bridges();
        vector<char> is_bridge(m, 0);
        for (auto bridge: bridges) is_bridge[bridge.id] = 1;

        int label = 0;
        vector<int> dfs_stack;
        fill(num.begin(), num.end(), -1);
        for (int root = 0; root < n; root++){
            if (num[root] != -1) continue;

            num[root] = label;
            dfs_stack.push_back(root);
            while (!dfs_stack.empty()){
                int u = dfs_stack.back();
                dfs_stack.pop_back();
                for (auto [v, id]: adj[u]){
                    if (num[v] == -1 && !is_bridge[id]){
                        num[v] = label;
                        dfs_stack.push_back(v);
                    }
                }
            }
            label++;
        }

        vector<Pair> bridge_tree;
        for (auto bridge: bridges){
            int u = bridge.u, v = bridge.v;
            bridge_tree.push_back(Pair(num[u], num[v]));
        }

        return bridge_tree;
    }
};

int main(){
    Graph graph(10);

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
    assert(graph.num == vector<int>({0, 0, 0, 0, 0, 1, 1, 1, 2, 3}));

    /// Closing the outer cycle 0-8-9-7 leaves no bridge, a second call relabels every node
    graph.add_edge(8, 0);
    graph.add_edge(9, 7);
    assert(graph.get_bridge_tree().empty());
    assert(graph.num == vector<int>(10, 0));

    Graph multi(4);
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
