/***
 *
 * Small Cycle Counting (triangles and 4-cycles)
 * Counts the 3-cycles and 4-cycles of a simple undirected graph, two cycles differ when their edge sets differ
 *
 * Complexity: O(n + m sqrt(m)) time, O(n + m) memory, for each count
 *
 * SmallCycles g(n); g.add_edge(u, v); with 0-based vertices
 *     the graph must be simple: no self loops (asserted) and no parallel edges (not checked, they inflate the counts)
 *     g.triangles() and g.four_cycles() return the counts as long long, both can be called any number of times
 *
 * Vertices are ranked by (degree, index) and every edge points from the lower rank to the higher one
 * A vertex has at most sqrt(2m) higher ranked neighbors, which bounds both scans
 * triangles: each triangle is found once, from its lowest ranked vertex along its two oriented edges
 * four_cycles: each 4-cycle is found once, from the vertex opposite its highest ranked vertex w,
 * by counting pairs of paths i - u - w with u -> w an oriented edge and w ranked above i
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SmallCycles{
    int n;
    vector<vector<int>> adj;

    SmallCycles(int n): n(n), adj(n) {}

    void add_edge(int u, int v){
        assert(u != v);
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    long long four_cycles() const{
        vector<vector<int>> up = oriented();
        vector<int> paths(n, 0);
        long long res = 0;

        for (int i = 0; i < n; i++){
            for (int u : adj[i]){
                for (int w : up[u]){
                    if (above(w, i)) res += paths[w]++;
                }
            }
            for (int u : adj[i]){
                for (int w : up[u]) paths[w] = 0;
            }
        }

        return res;
    }

    long long triangles() const{
        vector<vector<int>> up = oriented();
        vector<int> mark(n, -1);
        long long res = 0;

        for (int u = 0; u < n; u++){
            for (int v : up[u]) mark[v] = u;
            for (int v : up[u]){
                for (int w : up[v]) res += mark[w] == u;
            }
        }

        return res;
    }

private:
    bool above(int a, int b) const{
        return adj[a].size() != adj[b].size() ? adj[a].size() > adj[b].size() : a > b;
    }

    vector<vector<int>> oriented() const{
        vector<vector<int>> up(n);
        for (int u = 0; u < n; u++){
            for (int v : adj[u]){
                if (above(v, u)) up[u].push_back(v);
            }
        }
        return up;
    }
};

int main(){
    SmallCycles square(4);
    for (int i = 0; i < 4; i++) square.add_edge(i, (i + 1) % 4);
    assert(square.triangles() == 0);
    assert(square.four_cycles() == 1);

    /// K_4: C(4, 3) triangles and 3 C(4, 4) 4-cycles
    SmallCycles k4(4);
    for (int u = 0; u < 4; u++){
        for (int v = u + 1; v < 4; v++) k4.add_edge(u, v);
    }
    assert(k4.triangles() == 4);
    assert(k4.four_cycles() == 3);

    /// Wheel: hub 0 on the rim 1..5, a triangle per rim edge and a 4-cycle per rim path of three vertices
    SmallCycles wheel(6);
    for (int i = 1; i <= 5; i++) wheel.add_edge(0, i), wheel.add_edge(i, i % 5 + 1);
    assert(wheel.triangles() == 5);
    assert(wheel.four_cycles() == 5);
    return 0;
}
