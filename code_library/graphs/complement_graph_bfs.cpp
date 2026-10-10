/***
 *
 * Complement Graph BFS
 * BFS distances and connected components in the complement of a sparse undirected graph
 *
 * Complexity: O(n + m) per bfs() or components() call, O(n + m) memory, m = number of edges added
 *
 * ComplementGraph g(n); g.add_edge(u, v) for the edges of the ORIGINAL graph, vertices in [0, n)
 * The complement joins every pair u != v that has no original edge; self loops and duplicate edges are ignored
 * g.bfs(s): dist[v] = edges on a shortest complement path from s to v, -1 if unreachable
 * g.components(): comp[v] = complement component of v, numbered 0, 1, ... in order of their smallest vertex
 *
 * Every still unvisited vertex sits in one list: popping u scans the list, moves out the vertices
 * not adjacent to u and keeps the rest, which are original neighbours of u - so each scan costs
 * O(deg(u) + removed) and the complement's O(n^2) edges are never built
 * Iterative, so no recursion depth limit
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ComplementGraph{
    int n;
    vector<vector<int>> adj;

    ComplementGraph(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> bfs(int s) const{
        vector<int> dist(n, -1), unvisited, mark(n, -1), frontier;
        for (int v = 0; v < n; v++){
            if (v != s) unvisited.push_back(v);
        }

        dist[s] = 0;
        sweep(s, unvisited, mark, frontier, [&](int u, int v){ dist[v] = dist[u] + 1; });
        return dist;
    }

    vector<int> components() const{
        vector<int> comp(n, -1), unvisited(n), mark(n, -1), frontier;
        iota(unvisited.rbegin(), unvisited.rend(), 0);

        /// unvisited stays in decreasing order, so back() is the smallest unvisited vertex
        for (int count = 0; !unvisited.empty(); count++){
            int s = unvisited.back();
            unvisited.pop_back();
            comp[s] = count;
            sweep(s, unvisited, mark, frontier, [&](int, int v){ comp[v] = count; });
        }
        return comp;
    }

private:
    /// visit(u, v) fires once per newly reached v, with u its BFS parent
    template<typename F>
    void sweep(int s, vector<int>& unvisited, vector<int>& mark, vector<int>& frontier, F&& visit) const{
        frontier.assign(1, s);
        for (size_t i = 0; i < frontier.size() && !unvisited.empty(); i++){
            int u = frontier[i];
            for (int v : adj[u]) mark[v] = u;

            size_t kept = 0;
            for (int v : unvisited){
                if (mark[v] == u) unvisited[kept++] = v;
                else visit(u, v), frontier.push_back(v);
            }
            unvisited.resize(kept);
        }
    }
};

int main(){
    /***
     * The complement of the path 0 - 1 - 2 - 3 is the path 2 - 0 - 3 - 1
    ***/
    ComplementGraph path(4);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 3}}) path.add_edge(u, v);
    assert((path.bfs(2) == vector<int>{1, 3, 0, 2}));
    assert((path.components() == vector<int>{0, 0, 0, 0}));

    /***
     * The complement of the complete bipartite graph with parts {0, 3} and {1, 2, 4}
     * is the edge 0 - 3 plus the triangle 1 - 2 - 4
    ***/
    ComplementGraph bipartite(5);
    for (int u : {0, 3}){
        for (int v : {1, 2, 4}) bipartite.add_edge(u, v);
    }
    assert((bipartite.bfs(0) == vector<int>{0, -1, -1, 1, -1}));
    assert((bipartite.components() == vector<int>{0, 1, 1, 0, 1}));
    return 0;
}
