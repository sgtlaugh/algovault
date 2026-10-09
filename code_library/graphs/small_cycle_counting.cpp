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
    SmallCycles empty(0), single(1);
    assert(empty.triangles() == 0 && empty.four_cycles() == 0);
    assert(single.triangles() == 0 && single.four_cycles() == 0);

    SmallCycles square(4);
    for (int i = 0; i < 4; i++) square.add_edge(i, (i + 1) % 4);
    assert(square.triangles() == 0 && square.four_cycles() == 1);

    SmallCycles pentagon(5);
    for (int i = 0; i < 5; i++) pentagon.add_edge(i, (i + 1) % 5);
    assert(pentagon.triangles() == 0 && pentagon.four_cycles() == 0);

    /// C(n, 3) triangles and 3 C(n, 4) 4-cycles in K_n
    auto complete = [](int n){
        SmallCycles g(n);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) g.add_edge(u, v);
        }
        return g;
    };
    SmallCycles k3 = complete(3), k4 = complete(4), k5 = complete(5), k10 = complete(10);
    assert(k3.triangles() == 1 && k3.four_cycles() == 0);
    assert(k4.triangles() == 4 && k4.four_cycles() == 3);
    assert(k5.triangles() == 10 && k5.four_cycles() == 15);
    assert(k10.triangles() == 120 && k10.four_cycles() == 630);

    /// Wheel: hub 0 on the rim 1..5, a triangle per rim edge and a 4-cycle per rim path of three vertices
    SmallCycles wheel(6);
    for (int i = 1; i <= 5; i++) wheel.add_edge(0, i), wheel.add_edge(i, i % 5 + 1);
    assert(wheel.triangles() == 5 && wheel.four_cycles() == 5);

    /// K_{3,3}: C(3, 2)^2 4-cycles, no odd cycles
    SmallCycles k33(6);
    for (int a = 0; a < 3; a++){
        for (int b = 3; b < 6; b++) k33.add_edge(a, b);
    }
    assert(k33.triangles() == 0 && k33.four_cycles() == 9);

    /// Petersen graph has girth 5
    SmallCycles petersen(10);
    for (int i = 0; i < 5; i++){
        petersen.add_edge(i, (i + 1) % 5);
        petersen.add_edge(i, i + 5);
        petersen.add_edge(i + 5, (i + 2) % 5 + 5);
    }
    assert(petersen.triangles() == 0 && petersen.four_cycles() == 0);

    /// Hypercube Q_4: C(4, 2) 2^2 = 24 square faces, no triangles
    SmallCycles cube(16);
    for (int x = 0; x < 16; x++){
        for (int b = 0; b < 4; b++){
            if (x < (x ^ (1 << b))) cube.add_edge(x, x ^ (1 << b));
        }
    }
    assert(cube.triangles() == 0 && cube.four_cycles() == 24);

    /// K_{2,n} with n = 1e5: C(n, 2) = 4999950000 4-cycles overflows int
    int leaves = 100000;
    SmallCycles book(leaves + 2);
    for (int i = 2; i < leaves + 2; i++) book.add_edge(0, i), book.add_edge(1, i);
    assert(book.triangles() == 0 && book.four_cycles() == 4999950000LL);
    book.add_edge(0, 1);
    assert(book.triangles() == leaves && book.four_cycles() == 4999950000LL);

    return 0;
}
