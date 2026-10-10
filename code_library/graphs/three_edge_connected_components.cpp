/***
 *
 * 3-Edge-Connected Components (Tsin)
 * Groups the vertices of an undirected graph into classes where every pair is joined by 3 edge-disjoint paths
 * Equivalently, u and v share a class exactly when no removal of 2 or fewer edges disconnects them
 *
 * Complexity: O(n + m)
 *
 * ThreeEdgeConnected g(n); g.add_edge(u, v); auto groups = g.run();
 * Parallel edges and self loops are allowed, vertices are numbered from 0 to n - 1
 * g.count: number of components, g.comp[v]: component of v in [0, count)
 * Components are numbered by their smallest vertex, groups[c] lists the vertices of c in ascending order
 *
 * The DFS is iterative, so long paths cannot overflow the stack
 *
 * https://doi.org/10.1016/j.jda.2008.04.003 (Tsin, Yet another optimal algorithm for 3-edge-connectivity)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct ThreeEdgeConnected{
    int n, m = 0, count = 0;
    vector<vector<pair<int, int>>> adj; /// (neighbor, edge index)
    vector<int> comp;

    ThreeEdgeConnected(int n) : n(n), adj(n), comp(n, -1) {}

    void add_edge(int u, int v){
        adj[u].push_back({v, m});
        adj[v].push_back({u, m});
        m++;
    }

    vector<vector<int>> run(){
        in.assign(n, -1), out.assign(n, 0), low.assign(n, n), deg.assign(n, 0), path.assign(n, -1);
        nxt.resize(n);
        iota(nxt.begin(), nxt.end(), 0);
        vector<int> parent_edge(n, -1), it(n, 0), stack;
        int timer = 0;

        for (int s = 0; s < n; s++){
            if (in[s] != -1) continue;
            in[s] = timer++;
            stack.push_back(s);
            while (!stack.empty()){
                int v = stack.back();
                if (it[v] < (int)adj[v].size()){
                    auto [w, id] = adj[v][it[v]++];
                    if (w == v || id == parent_edge[v]) continue;
                    if (in[w] == -1){
                        in[w] = timer++, parent_edge[w] = id;
                        stack.push_back(w);
                    }
                    else if (in[w] < in[v]){
                        deg[v]++;
                        low[v] = min(low[v], in[w]);
                    }
                    else close_cycle(v, w);
                    continue;
                }

                out[v] = timer;
                stack.pop_back();
                if (!stack.empty()) merge_child(stack.back(), v);
            }
        }

        count = 0;
        fill(comp.begin(), comp.end(), -1);
        for (int s = 0; s < n; s++){
            if (comp[s] != -1) continue;
            for (int v = s; comp[v] == -1; v = nxt[v]) comp[v] = count;
            count++;
        }

        vector<vector<int>> groups(count);
        for (int v = 0; v < n; v++) groups[comp[v]].push_back(v);
        return groups;
    }

private:
    /// path[v]: the chain of descendants that will join v's class once a back edge closes over them
    /// nxt: classes kept as cyclic lists, merged in O(1) by swapping successors
    vector<int> in, out, low, deg, path, nxt;

    void absorb(int v, int w){
        swap(nxt[v], nxt[w]);
        deg[v] += deg[w];
    }

    /// Back edge (w, v) seen from its upper end v: every chain vertex whose subtree holds w joins v
    void close_cycle(int v, int w){
        deg[v]--;
        int u = path[v];
        while (u != -1 && in[u] <= in[w] && in[w] < out[u]){
            absorb(v, u);
            u = path[u];
        }
        path[v] = u;
    }

    void merge_child(int v, int w){
        if (path[w] == -1 && deg[w] <= 1){
            deg[v] += deg[w];
            low[v] = min(low[v], low[w]);
            return;
        }
        if (deg[w] == 0) w = path[w];

        /// On equal low either chain may stay as path[v], the other is absorbed into v either way
        if (low[w] < low[v]){
            low[v] = low[w];
            swap(w, path[v]);
        }
        for (; w != -1; w = path[w]) absorb(v, w);
    }
};

int main(){
    /***
     *   0 ----- 1
     *   | \     |
     *   |   \   |
     *   |     \ |
     *   3 ----- 2
    ***/
    ThreeEdgeConnected g(4);
    g.add_edge(0, 2);
    g.add_edge(0, 1);
    g.add_edge(3, 0);
    g.add_edge(2, 1);
    g.add_edge(2, 3);
    assert((g.run() == vector<vector<int>>{{0, 2}, {1}, {3}}));  /// 1 and 3 have only 2 edges each
    assert(g.count == 3);
    assert((g.comp == vector<int>{0, 1, 0, 2}));

    /// Removing 2 edges splits a cycle, a doubled cycle needs 4
    ThreeEdgeConnected cycle(5), doubled(5);
    for (int i = 0; i < 5; i++){
        cycle.add_edge(i, (i + 1) % 5);
        doubled.add_edge(i, (i + 1) % 5), doubled.add_edge(i, (i + 1) % 5);
    }
    assert(cycle.run().size() == 5);
    assert(doubled.run().size() == 1);
    return 0;
}
