/***
 *
 * Euler Path (Hierholzer)
 * Finds a walk that uses every edge exactly once, in a directed or undirected multigraph
 *
 * Complexity: O(n + m)
 *
 * EulerDirected g(n) or EulerUndirected g(n), vertices numbered from 0 to n-1
 * g.add_edge(u, v) returns the edge id, ids are 0, 1, 2, ... in call order
 * g.solve(start) returns the m + 1 vertices of an Euler path beginning at start, or an empty vector if none exists
 * start must be -1 or a vertex in [0, n - 1]
 * After a successful solve, g.walk_edges holds the m edge ids in order: edge walk_edges[i] joins path[i] to path[i + 1]
 * start = -1 picks a valid start: the vertex with out - in = 1 (directed) or the first odd degree vertex (undirected),
 * otherwise the first vertex with an outgoing edge
 * The path is a cycle iff path.front() == path.back()
 * With no edges the path is the single vertex start (0 for start = -1), and it is empty when n = 0
 *
 * Multi-edges, self loops and isolated vertices are allowed: only the vertices that have edges need to be connected
 * Iterative, so long paths cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<bool directed>
struct EulerPath{
    int n, m = 0;
    vector<vector<pair<int, int>>> adj;
    vector<int> balance, walk_edges; /// balance: out - in when directed, the degree when undirected

    EulerPath(int n) : n(n), adj(n), balance(n, 0) {}

    int add_edge(int u, int v){
        adj[u].push_back({v, m});
        if (directed) balance[u]++, balance[v]--;
        else adj[v].push_back({u, m}), balance[u]++, balance[v]++;
        return m++;
    }

    vector<int> solve(int start = -1){
        walk_edges.clear();
        if (n == 0) return {};

        int unbalanced = 0, source = -1;
        for (int u = 0; u < n; u++){
            int b = excess(u);
            if (b == 0) continue;
            if (b != 1 && b != -1) return {};
            unbalanced++;
            if (b == 1 && source == -1) source = u;
        }
        if (unbalanced > 2) return {};

        if (start == -1){
            start = unbalanced ? source : 0;
            for (int u = 0; u < n && !unbalanced; u++){
                if (!adj[u].empty()){
                    start = u;
                    break;
                }
            }
        }
        if (unbalanced ? excess(start) != 1 : m > 0 && adj[start].empty()) return {};

        vector<int> it(n, 0), path;
        vector<char> used(m, 0);
        vector<pair<int, int>> stack = {{start, -1}};
        while (!stack.empty()){
            auto [u, e] = stack.back();
            if (it[u] == (int)adj[u].size()){
                path.push_back(u);
                if (e != -1) walk_edges.push_back(e);
                stack.pop_back();
                continue;
            }
            auto [v, id] = adj[u][it[u]++];
            if (!used[id]) used[id] = 1, stack.push_back({v, id});
        }

        if ((int)path.size() != m + 1){
            walk_edges.clear();
            return {};
        }
        reverse(path.begin(), path.end());
        reverse(walk_edges.begin(), walk_edges.end());
        return path;
    }

private:
    int excess(int u) const{
        return directed ? balance[u] : balance[u] % 2;
    }
};

using EulerDirected = EulerPath<true>;
using EulerUndirected = EulerPath<false>;

int main(){
    /***
     * 0 -> 1 -> 2 -> 0 and 0 -> 3: vertex 0 has out - in = 1, and taking 0 -> 3 first strands the cycle
    ***/
    EulerDirected d(4);
    d.add_edge(0, 1), d.add_edge(1, 2), d.add_edge(2, 0);
    assert(d.add_edge(0, 3) == 3);                       /// ids follow the call order
    assert((d.solve() == vector<int>{0, 1, 2, 0, 3}));
    assert((d.walk_edges == vector<int>{0, 1, 2, 3}));
    assert(d.solve(3).empty());                          /// every Euler path here leaves from 0

    /***
     * Path 0 - 1 - 2 with a self loop at 1: the odd vertices are 0 and 2, vertex 1 has degree 4
    ***/
    EulerUndirected u(3);
    u.add_edge(0, 1), u.add_edge(1, 1), u.add_edge(1, 2);
    assert((u.solve(2) == vector<int>{2, 1, 1, 0}));
    assert((u.walk_edges == vector<int>{2, 1, 0}));
    assert(u.solve(1).empty());                          /// must start at an odd vertex

    EulerUndirected star(4);
    star.add_edge(0, 1), star.add_edge(0, 2), star.add_edge(0, 3);
    assert(star.solve().empty());                        /// four odd vertices
    return 0;
}
