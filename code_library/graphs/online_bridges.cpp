/***
 *
 * Online Bridges
 * Maintains the bridge count and the 2-edge-connected components of an undirected graph under edge insertions
 *
 * Complexity: O((n log n + m) α(n)) total for m insertions, O(n) memory
 *
 * OnlineBridges g(n): an empty graph on vertices 0 to n - 1
 * g.add_edge(u, v): adds an edge, parallel edges and self loops allowed
 * g.component(u): representative of the 2-edge-connected component of u, changes as components merge
 * g.connected(u, v): u and v are in the same connected component
 * g.two_edge_connected(u, v): u and v are in the same 2-edge-connected component
 * g.bridges: number of bridges in the current graph
 *
 * Notes:
 *   - An existing edge (u, v) is a bridge exactly when !g.two_edge_connected(u, v)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct OnlineBridges{
    int bridges = 0, iteration = 0;
    vector<int> two_edge_link, rnk;
    vector<int> par;  /// at a representative: a vertex of its parent component in the bridge tree, -1 at roots
    vector<int> conn_link, conn_size;
    vector<int> last_visit;

    OnlineBridges(int n) : two_edge_link(n), rnk(n, 0), par(n, -1), conn_link(n), conn_size(n, 1), last_visit(n, 0){
        iota(two_edge_link.begin(), two_edge_link.end(), 0);
        iota(conn_link.begin(), conn_link.end(), 0);
    }

    void add_edge(int u, int v){
        int a = component(u), b = component(v);

        if (a == b) return;

        int ca = find(conn_link, a), cb = find(conn_link, b);
        if (ca == cb){
            merge_path(a, b);
            return;
        }

        bridges++;
        if (conn_size[ca] > conn_size[cb]) swap(a, b), swap(ca, cb);
        make_root(a);
        par[a] = b;
        conn_link[ca] = cb;
        conn_size[cb] += conn_size[ca];
    }

    int component(int u){
        return find(two_edge_link, u);
    }

    bool connected(int u, int v){
        return find(conn_link, u) == find(conn_link, v);
    }

    bool two_edge_connected(int u, int v){
        return component(u) == component(v);
    }

    int find(vector<int>& link, int u){
        int root = u;

        while (link[root] != root) root = link[root];

        while (link[u] != root){
            int next = link[u];
            link[u] = root;
            u = next;
        }

        return root;
    }

    /// Walks up from both ends in turns, the first component seen twice is the LCA
    void merge_path(int a, int b){
        vector<int> path_a, path_b;
        int lca = -1;

        iteration++;
        auto step = [&](int& x, vector<int>& path){
            if (x == -1 || lca != -1) return;
            x = component(x);
            path.push_back(x);
            if (last_visit[x] == iteration) lca = x;
            else last_visit[x] = iteration, x = par[x];
        };
        while (lca == -1) step(a, path_a), step(b, path_b);

        int lca_par = par[lca];
        for (auto path : {&path_a, &path_b}){
            for (int x : *path){
                if (x == lca) break;
                unite(x, lca);
                bridges--;
            }
        }
        par[component(lca)] = lca_par;
    }

    void unite(int x, int y){
        x = component(x), y = component(y);
        if (rnk[x] < rnk[y]) swap(x, y);
        two_edge_link[y] = x;
        if (rnk[x] == rnk[y]) rnk[x]++;
    }

    void make_root(int u){
        int child = -1;
        while (u != -1){
            int up = par[u] == -1 ? -1 : component(par[u]);
            par[u] = child;
            child = u;
            u = up;
        }
    }
};

int main(){
    /// 0 - 1 - 2 - 0 is a triangle, 3 - 4 - 5 a path that 2 - 3 attaches and 5 - 3 closes into a second triangle
    OnlineBridges g(6);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    assert(g.bridges == 2);
    g.add_edge(2, 0);  /// closes the triangle, both bridges disappear
    assert(g.bridges == 0);
    assert(g.two_edge_connected(0, 2));

    g.add_edge(3, 4);
    g.add_edge(4, 5);
    assert(!g.connected(0, 3));
    g.add_edge(2, 3);
    assert(g.connected(0, 3));
    assert(g.bridges == 3);

    g.add_edge(5, 3);
    assert(g.bridges == 1);  /// only 2 - 3 is left
    assert(!g.two_edge_connected(2, 3));
    assert(g.two_edge_connected(3, 5));
    return 0;
}
