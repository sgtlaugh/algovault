/***
 *
 * Minimum Spanning Tree (Kruskal, Boruvka, dense Prim)
 * Minimum spanning forest of an undirected weighted graph
 *
 * Complexity: Kruskal O(m log m), Boruvka O((n + m) log n), dense_prim O(n^2)
 *
 * Kruskal mst(n): an empty graph on vertices 0 to n - 1, Boruvka mst(n) has the same API
 * mst.add_edge(u, v, w): adds an edge of long long weight w, negative allowed
 * mst.solve(): {total weight, indices of the chosen edges in add_edge order, connected}
 * dense_prim(n, w): {total weight, parent, connected} for edge weights w(u, v) of type int or long long
 *
 * Notes:
 *   - connected is true when the forest is one tree spanning all n vertices, else it has one tree per component
 *   - Kruskal lists chosen edges by increasing weight, Boruvka in the order it picks them
 *   - The total weight of the forest must fit in long long
 *   - For implicit graphs like xor MST, swap Boruvka's edge scan for a query of a lightest edge leaving each component
 *   - dense_prim: w(u, v) == numeric_limits<T>::max() means no edge for the type T that w returns
 *   - dense_prim: w is called once per pair {u, v} in either order and never with u == v, so it must be symmetric
 *   - dense_prim: parent[v] is v's neighbour in the forest, -1 at one unspecified root per tree
 *   - dense_prim is for complete graphs where Kruskal's sort is the bottleneck, w computed on the fly needs no matrix
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct MSTEdge{
    int u, v;
    long long w;
};

struct MSTResult{
    long long weight;
    vector<int> chosen;
    bool connected;
};

struct DenseMSTResult{
    long long weight;
    vector<int> parent;
    bool connected;
};

struct MSTForest{
    vector<int> parent, rank;

    MSTForest(int n) : parent(n), rank(n, 0){
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    }

    bool connect(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (rank[a] < rank[b]) swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
        return true;
    }
};

struct Kruskal{
    int n;
    vector<MSTEdge> edges;

    Kruskal(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    MSTResult solve() const{
        vector<int> order(edges.size());

        iota(order.begin(), order.end(), 0);
        stable_sort(order.begin(), order.end(), [&](int a, int b){ return edges[a].w < edges[b].w; });

        MSTForest dsu(n);
        MSTResult res{0, {}, false};
        for (int i : order){
            if (!dsu.connect(edges[i].u, edges[i].v)) continue;
            res.weight += edges[i].w;
            res.chosen.push_back(i);
        }

        res.connected = (int)res.chosen.size() == max(0, n - 1);
        return res;
    }
};

struct Boruvka{
    int n;
    vector<MSTEdge> edges;

    Boruvka(int n) : n(n) {}

    void add_edge(int u, int v, long long w){
        edges.push_back({u, v, w});
    }

    MSTResult solve() const{
        MSTForest dsu(n);
        MSTResult res{0, {}, false};
        vector<int> best(n);

        for (bool merged = true; merged; ){
            merged = false;
            fill(best.begin(), best.end(), -1);
            for (int i = 0; i < (int)edges.size(); i++){
                int a = dsu.find(edges[i].u), b = dsu.find(edges[i].v);
                if (a == b) continue;
                if (best[a] == -1 || edges[i].w < edges[best[a]].w) best[a] = i;
                if (best[b] == -1 || edges[i].w < edges[best[b]].w) best[b] = i;
            }

            /// both endpoints can pick an edge and implicit queries break ties arbitrarily: connect() drops extra picks
            for (int c = 0; c < n; c++){
                int i = best[c];
                if (i == -1 || !dsu.connect(edges[i].u, edges[i].v)) continue;
                res.weight += edges[i].w;
                res.chosen.push_back(i);
                merged = true;
            }
        }

        res.connected = (int)res.chosen.size() == max(0, n - 1);
        return res;
    }
};

template<typename Weight>
DenseMSTResult dense_prim(int n, Weight w){
    using T = decltype(w(0, 0));
    const T NONE = numeric_limits<T>::max();
    int roots = 0;
    vector<T> best(n, NONE);
    vector<char> used(n, 0);
    DenseMSTResult res{0, vector<int>(n, -1), false};

    for (int it = 0, u = 0; it < n; it++){
        if (best[u] == NONE) roots++;
        else res.weight += best[u];
        used[u] = 1;

        /// Picking the next vertex in the same pass as the update skips a second scan, 1.5x faster at n = 5000
        int next = -1;
        for (int v = 0; v < n; v++){
            if (used[v]) continue;
            T cost = w(u, v);
            if (cost < best[v]) best[v] = cost, res.parent[v] = u;
            if (next == -1 || best[v] < best[next]) next = v;
        }
        u = next;
    }

    res.connected = roots <= 1;
    return res;
}

int main(){
    /***
     *   0 --1-- 1
     *   |  \    |
     *   4   3   2
     *   |    \  |
     *   3 --5-- 2
    ***/
    vector<MSTEdge> edges = {{0, 1, 1}, {1, 2, 2}, {0, 2, 3}, {0, 3, 4}, {3, 2, 5}};
    Kruskal kruskal(4);
    Boruvka boruvka(4);
    for (auto [u, v, w] : edges) kruskal.add_edge(u, v, w), boruvka.add_edge(u, v, w);

    MSTResult r = kruskal.solve();
    assert(r.weight == 7);
    assert(r.connected);
    assert((r.chosen == vector<int>{0, 1, 3}));     /// edge indices by increasing weight
    assert(boruvka.solve().weight == 7);

    /// Two components give a spanning forest, negative weights are fine
    Kruskal forest(4);
    forest.add_edge(0, 1, -4), forest.add_edge(2, 3, 6);
    assert(forest.solve().weight == 2);
    assert(!forest.solve().connected);

    const long long N = LLONG_MAX;   /// no edge
    vector<vector<long long>> w = {{N, 1, 3, 4}, {1, N, 2, N}, {3, 2, N, 5}, {4, N, 5, N}};
    DenseMSTResult d = dense_prim(4, [&](int u, int v){ return w[u][v]; });
    assert(d.weight == 7);
    assert((d.parent == vector<int>{-1, 0, 1, 0}));  /// the same tree as above, rooted at 0

    /// Complete graph on points with Manhattan distances computed on the fly
    vector<pair<int, int>> points = {{0, 0}, {3, 0}, {3, 4}, {0, 4}, {1, 1}};
    auto manhattan = [&](int u, int v){
        return abs(points[u].first - points[v].first) + abs(points[u].second - points[v].second);
    };
    assert(dense_prim(5, manhattan).weight == 12);
    return 0;
}
