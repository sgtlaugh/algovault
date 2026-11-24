/***
 *
 * Find Cycle
 * Any one cycle of a directed or undirected graph, or the fact that it has none
 *
 * Complexity: O(n + m)
 *
 * FindCycle g(n, FindCycle::DIRECTED): empty graph on vertices 0 to n - 1, or UNDIRECTED
 * g.add_edge(u, v): edge u -> v, or u - v, its id counts up from 0 in add_edge order
 * g.solve(): {exists, vertices, edges}, one simple cycle, both lists empty when the graph has none
 *
 * Notes:
 *   - vertices[i] -> vertices[i + 1] uses edges[i], and the last edge returns to vertices[0]
 *   - A self loop is a cycle of one; undirected cycles never reuse an edge, so two parallel edges are a cycle of two
 *   - No cycle means a DAG when directed, a forest when undirected
 *   - Iterative, so 5e5 vertices on one path cannot overflow the stack
 *   - For the lightest cycle use minimum_weight_cycle.cpp, which costs O(n m)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FindCycle{
    enum Kind{DIRECTED, UNDIRECTED};

    struct Cycle{
        bool exists;
        vector<int> vertices, edges;
    };

    int n;
    Kind kind;
    vector<vector<pair<int, int>>> adj;
    int m = 0;

    FindCycle(int n, Kind kind) : n(n), kind(kind), adj(n) {}

    void add_edge(int u, int v){
        assert(0 <= u && u < n && 0 <= v && v < n);
        adj[u].push_back({v, m});
        if (kind == UNDIRECTED && u != v) adj[v].push_back({u, m});
        m++;
    }

    Cycle solve() const{
        /// state 1 marks the vertices on the DFS stack: an edge into one of them closes a cycle
        vector<int> state(n, 0), from(n, -1), from_edge(n, -1), it(n, 0), stack;
        for (int s = 0; s < n; s++){
            if (state[s]) continue;

            state[s] = 1, stack.push_back(s);
            while (!stack.empty()){
                int u = stack.back();
                if (it[u] == (int)adj[u].size()){
                    state[u] = 2, stack.pop_back();
                    continue;
                }

                auto [v, id] = adj[u][it[u]++];
                if (kind == UNDIRECTED && id == from_edge[u]) continue;
                if (state[v] == 1) return unwind(u, v, id, from, from_edge);
                if (state[v] == 0){
                    state[v] = 1, from[v] = u, from_edge[v] = id;
                    stack.push_back(v);
                }
            }
        }

        return {false, {}, {}};
    }

    /// The tree path v -> ... -> u, then the closing edge u -> v
    static Cycle unwind(int u, int v, int id, const vector<int>& from, const vector<int>& from_edge){
        Cycle res{true, {}, {}};
        for (int x = u; x != v; x = from[x]){
            res.vertices.push_back(x);
            res.edges.push_back(from_edge[x]);
        }
        res.vertices.push_back(v);
        reverse(res.vertices.begin(), res.vertices.end());
        reverse(res.edges.begin(), res.edges.end());
        res.edges.push_back(id);

        return res;
    }
};

int main(){
    /// Directed: 0 -> 1 -> 2 -> 0 is the only cycle, 3 hangs off it; the DFS from 0 closes it at edge 2
    FindCycle d(4, FindCycle::DIRECTED);
    d.add_edge(0, 1);
    d.add_edge(1, 2);
    d.add_edge(2, 0);
    d.add_edge(2, 3);
    FindCycle::Cycle c = d.solve();
    assert(c.exists);
    assert((c.vertices == vector<int>{0, 1, 2}));
    assert((c.edges == vector<int>{0, 1, 2}));

    /// A DAG has no cycle
    FindCycle dag(3, FindCycle::DIRECTED);
    dag.add_edge(0, 1);
    dag.add_edge(0, 2);
    dag.add_edge(1, 2);
    assert(!dag.solve().exists);

    /// Undirected: the parallel edges 0 and 1 between vertices 0 and 1 form a cycle of two
    FindCycle u(3, FindCycle::UNDIRECTED);
    u.add_edge(0, 1);
    u.add_edge(1, 0);
    u.add_edge(1, 2);
    FindCycle::Cycle two = u.solve();
    assert((two.vertices == vector<int>{0, 1}));
    assert((two.edges == vector<int>{0, 1}));

    /// A self loop is a cycle of one
    FindCycle loop(2, FindCycle::UNDIRECTED);
    loop.add_edge(0, 1);
    loop.add_edge(1, 1);
    assert((loop.solve().vertices == vector<int>{1}));

    return 0;
}
