/***
 *
 * Chordal Graph
 * Tests if every cycle of 4 or more vertices has a chord, proven by a perfect elimination order or an induced cycle
 *
 * Complexity: O(n + m)
 *
 * ChordalGraph g(n): an empty undirected graph on vertices 0 to n - 1
 * g.add_edge(u, v): adds the edge u - v
 * g.run(): true if the graph is chordal, then fills order, otherwise fills cycle
 * g.max_clique(): the vertices of one maximum clique, after run() returned true
 * g.max_independent_set(): the vertices of one maximum independent set, after run() returned true
 * g.order: all vertices, each forming a clique with its neighbors later in order, after run() returned true
 * g.cycle: the vertices of a chordless cycle of length at least 4 in cycle order, empty when chordal
 *
 * Notes:
 *   - The graph must be simple: u != v and no edge added twice
 *   - run() can be called again after adding edges
 *   - Chordal graphs are perfect, so the chromatic number is max_clique().size()
 *   - Maximum cardinality search with the certificate of Tarjan and Yannakakis, SIAM J. Comput. 13 (1984)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ChordalGraph{
    int n;
    vector<vector<int>> adj;
    vector<int> order, cycle;

    ChordalGraph(int n) : n(n), adj(n) {}

    void add_edge(int u, int v){
        assert(u != v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    bool run(){
        order = maximum_cardinality_search();
        cycle.clear();

        vector<int> pos(n), parent(n, -1);
        for (int i = 0; i < n; i++) pos[order[i]] = i;
        for (int v = 0; v < n; v++){
            for (int w : adj[v]){
                if (pos[w] > pos[v] && (parent[v] == -1 || pos[w] < pos[parent[v]])) parent[v] = w;
            }
        }

        /// v's later neighbors form a clique iff those besides parent[v] are adjacent to parent[v]
        vector<int> mark(n, -1);
        for (int p = 0; p < n; p++){
            for (int x : adj[p]) mark[x] = p;
            for (int v : adj[p]){
                if (parent[v] != p) continue;
                for (int w : adj[v]){
                    if (pos[w] > pos[p] && mark[w] != p){
                        cycle = induced_cycle(v, p, w);
                        return false;
                    }
                }
            }
        }

        return true;
    }

    vector<int> max_clique() const{
        assert((int)order.size() == n && cycle.empty());

        vector<int> pos(n);
        for (int i = 0; i < n; i++) pos[order[i]] = i;

        int best = -1, best_size = 0;
        for (int v = 0; v < n; v++){
            int size = 1;
            for (int w : adj[v]) size += pos[w] > pos[v];
            if (size > best_size) best = v, best_size = size;
        }

        vector<int> clique;
        if (best == -1) return clique;
        clique.push_back(best);
        for (int w : adj[best]){
            if (pos[w] > pos[best]) clique.push_back(w);
        }

        return clique;
    }

    vector<int> max_independent_set() const{
        assert((int)order.size() == n && cycle.empty());

        /// Taking each vertex that is simplicial in what remains is optimal (Gavril 1972)
        vector<int> res;
        vector<char> blocked(n, 0);
        for (int v : order){
            if (blocked[v]) continue;
            res.push_back(v);
            for (int w : adj[v]) blocked[w] = 1;
        }

        return res;
    }

    /// The reverse of the visit order, so the first vertex in order is the last one visited
    vector<int> maximum_cardinality_search() const{
        vector<int> weight(n, 0), visit;
        vector<vector<int>> bucket(n + 1);
        vector<char> visited(n, 0);
        for (int v = 0; v < n; v++) bucket[0].push_back(v);

        int best = 0;
        while ((int)visit.size() < n){
            while (bucket[best].empty()) best--;
            int v = bucket[best].back();
            bucket[best].pop_back();
            /// No bucket above best holds entries, so a lower stale copy surfaces only after its vertex was visited
            if (visited[v]) continue;

            visited[v] = 1;
            visit.push_back(v);
            for (int w : adj[v]){
                if (visited[w]) continue;
                bucket[++weight[w]].push_back(w);
                best = max(best, weight[w]);
            }
        }

        reverse(visit.begin(), visit.end());
        return visit;
    }

    /// Shortest p - w path avoiding v's other neighbors, closed by v; MCS guarantees the path exists
    vector<int> induced_cycle(int v, int p, int w) const{
        vector<int> from(n, -1);
        vector<char> blocked(n, 0);
        blocked[v] = 1;
        for (int x : adj[v]) blocked[x] = 1;
        blocked[w] = 0;

        vector<int> queue = {p};
        from[p] = p;
        for (int i = 0; i < (int)queue.size() && from[w] == -1; i++){
            for (int x : adj[queue[i]]){
                if (!blocked[x] && from[x] == -1) from[x] = queue[i], queue.push_back(x);
            }
        }
        assert(from[w] != -1);

        vector<int> res = {v};
        for (int x = w; x != p; x = from[x]) res.push_back(x);
        res.push_back(p);

        return res;
    }
};

int main(){
    /// Two triangles 0-1-2 and 1-2-3 sharing edge 1-2, plus the pendant 4 on 3
    ChordalGraph diamond(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}, {3, 4}}) diamond.add_edge(u, v);
    assert(diamond.run());
    assert(diamond.cycle.empty());
    assert(diamond.max_clique().size() == 3);
    assert(diamond.max_independent_set().size() == 2);  /// at most one from triangle 0-1-2 and one from edge 3-4

    /// The square 0-1-2-3 has no chord, so it is its own certificate
    ChordalGraph square(4);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 3}, {3, 0}}) square.add_edge(u, v);
    assert(!square.run());
    assert(square.cycle.size() == 4);

    square.add_edge(0, 2);
    assert(square.run());
    assert(square.max_clique().size() == 3);
    return 0;
}
