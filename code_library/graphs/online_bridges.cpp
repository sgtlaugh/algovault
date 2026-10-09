/***
 *
 * Online Bridges
 * Maintains the bridge count and the 2-edge-connected components of an undirected graph under edge insertions
 *
 * Complexity: O((n log n + m) α(n)) total for m insertions, O(n) memory
 *
 * OnlineBridges g(n); g.add_edge(u, v);   vertices 0 to n - 1, parallel edges and self loops allowed
 * g.bridges: number of bridges in the current graph
 * g.component(u): representative of the 2-edge-connected component of u, changes as components merge
 * g.two_edge_connected(u, v): u and v are in the same 2-edge-connected component
 * g.connected(u, v): u and v are in the same connected component
 * An existing edge (u, v) is a bridge exactly when !g.two_edge_connected(u, v)
 *
 * Each connected component is kept as a rooted tree of its 2-edge-connected components (the bridge tree)
 * An edge between two trees re-roots the smaller one and hangs it below the other, O(n log n) total
 * An edge inside one tree closes a cycle, every component on it merges into their LCA and its bridges disappear
 * cp-algorithms merges without union by rank for O(n log n + m log n), the rank here keeps finds at α(n)
 *
 * Everything is iterative, so long paths cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct OnlineBridges{
    int bridges = 0, iteration = 0;
    vector<int> dsu, rnk;  /// 2-edge-connected components
    vector<int> par;  /// at a component representative: any vertex of its parent component in the bridge tree, -1 at a root
    vector<int> cc, cc_size;  /// connected components
    vector<int> last_visit;

    OnlineBridges(int n) : dsu(n), rnk(n, 0), par(n, -1), cc(n), cc_size(n, 1), last_visit(n, 0){
        iota(dsu.begin(), dsu.end(), 0);
        iota(cc.begin(), cc.end(), 0);
    }

    void add_edge(int u, int v){
        int a = component(u), b = component(v);
        if (a == b) return;

        int ca = find_cc(a), cb = find_cc(b);
        if (ca == cb){
            merge_path(a, b);
            return;
        }

        bridges++;
        if (cc_size[ca] > cc_size[cb]) swap(a, b), swap(ca, cb);
        make_root(a);
        par[a] = b;
        cc[ca] = cb;
        cc_size[cb] += cc_size[ca];
    }

    int component(int u){
        int root = u;
        while (dsu[root] != root) root = dsu[root];

        while (dsu[u] != root){
            int next = dsu[u];
            dsu[u] = root;
            u = next;
        }
        return root;
    }

    bool connected(int u, int v){
        return find_cc(u) == find_cc(v);
    }

    bool two_edge_connected(int u, int v){
        return component(u) == component(v);
    }

    int find_cc(int u){
        int root = u;
        while (cc[root] != root) root = cc[root];

        while (cc[u] != root){
            int next = cc[u];
            cc[u] = root;
            u = next;
        }
        return root;
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
        dsu[y] = x;
        if (rnk[x] == rnk[y]) rnk[x]++;
    }
};

int main(){
    /***
     * 0 - 1 - 2 - 0 is a triangle, 2 - 3 and 4 - 5 hang off it, 5 - 3 closes the triangle 3 - 4 - 5, 0 - 5 closes everything
    ***/
    OnlineBridges g(6);
    vector<pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {2, 3}, {4, 5}, {5, 3}, {0, 5}};
    vector<int> expected = {1, 2, 0, 1, 2, 3, 1, 0};
    for (int i = 0; i < (int)edges.size(); i++){
        g.add_edge(edges[i].first, edges[i].second);
        assert(g.bridges == expected[i]);
        if (i == 3) assert(!g.connected(0, 3) && g.connected(3, 4) && g.two_edge_connected(0, 2) && !g.two_edge_connected(3, 4));
        if (i == 6) assert(g.connected(0, 5) && !g.two_edge_connected(2, 3) && g.two_edge_connected(3, 5));
    }
    for (int v = 1; v < 6; v++) assert(g.component(v) == g.component(0));

    /***
     * Chain 0 - 1 - 2 and edge 3 - 4 join through 4 - 1, then 3 - 0 closes the cycle 0 - 1 - 4 - 3, leaving 1 - 2 the only bridge
    ***/
    OnlineBridges rerooted(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {3, 4}, {4, 1}}) rerooted.add_edge(u, v);
    assert(rerooted.bridges == 4);
    rerooted.add_edge(3, 0);
    assert(rerooted.bridges == 1);
    assert(rerooted.two_edge_connected(0, 4) && rerooted.two_edge_connected(1, 3) && !rerooted.two_edge_connected(1, 2));

    OnlineBridges parallel(3);
    parallel.add_edge(0, 1);
    parallel.add_edge(1, 2);
    assert(parallel.bridges == 2);
    parallel.add_edge(1, 0);
    assert(parallel.bridges == 1 && parallel.two_edge_connected(0, 1) && !parallel.two_edge_connected(1, 2));
    parallel.add_edge(2, 2);
    assert(parallel.bridges == 1);

    OnlineBridges single(1);
    single.add_edge(0, 0);
    assert(single.bridges == 0 && single.component(0) == 0);

    return 0;
}
