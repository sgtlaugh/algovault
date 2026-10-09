/***
 *
 * Minimum Spanning Tree (Kruskal, Boruvka, dense Prim)
 * Minimum spanning forest of an undirected weighted graph
 *
 * Complexity: Kruskal O(m log m), Boruvka O((n + m) log n), dense_prim O(n^2)
 *
 * Kruskal and Boruvka share one API: 0-based vertices, long long weights (negative allowed)
 *     Kruskal mst(n); mst.add_edge(u, v, w); MSTResult r = mst.solve();
 *     the total weight of the forest must fit in long long
 * solve() returns {total weight, indices of the chosen edges in add_edge order, connected}
 *     connected is true when the forest is one tree spanning all n vertices
 *     for a disconnected graph the result is a minimum spanning forest, one tree per component
 *     Kruskal lists chosen edges by increasing weight, Boruvka in the order it picks them
 *
 * Boruvka is the template for implicit graphs (xor MST, Manhattan MST, complete graphs given by a formula)
 *     each round only needs the cheapest edge leaving every component: replace the edge scan in solve()
 *     with a per-component query, any lightest edge per component works since ties are resolved by the union-find
 *     there are at most log2(n) + 1 rounds since every round at least halves the number of components
 *     that still have an outgoing edge, and the last round merges nothing
 *
 * dense_prim(w): w is an n x n symmetric matrix of int or long long, w[u][v] == numeric_limits<T>::max() means no edge
 *     for complete graphs where m ~ n^2 makes Kruskal's sort the bottleneck, n ~ 5000 with an int matrix (100 MB)
 *     returns {total weight, parent, connected}: parent[v] is v's neighbour in the forest, -1 for one root per tree
 *     the total weight of the forest must fit in long long
 *     the diagonal is ignored
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

    MSTForest(int n) : parent(n), rank(n, 0) {
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int x){
        while (x != parent[x]) x = parent[x] = parent[parent[x]];
        return x;
    }

    bool unite(int a, int b){
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
            if (!dsu.unite(edges[i].u, edges[i].v)) continue;
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

            /// both endpoints of an edge can pick it, and implicit-graph queries may break ties arbitrarily: unite() drops the redundant picks
            for (int c = 0; c < n; c++){
                int i = best[c];
                if (i == -1 || !dsu.unite(edges[i].u, edges[i].v)) continue;
                res.weight += edges[i].w;
                res.chosen.push_back(i);
                merged = true;
            }
        }

        res.connected = (int)res.chosen.size() == max(0, n - 1);
        return res;
    }
};

template<typename T>
DenseMSTResult dense_prim(const vector<vector<T>>& w){
    const T NONE = numeric_limits<T>::max();
    int n = w.size(), roots = 0;
    vector<T> best(n, NONE);
    vector<char> used(n, 0);
    DenseMSTResult res{0, vector<int>(n, -1), false};

    for (int it = 0; it < n; it++){
        int u = -1;
        for (int v = 0; v < n; v++){
            if (!used[v] && (u == -1 || best[v] < best[u])) u = v;
        }

        if (best[u] == NONE) roots++;
        else res.weight += best[u];
        used[u] = 1;

        for (int v = 0; v < n; v++){
            if (!used[v] && w[u][v] < best[v]) best[v] = w[u][v], res.parent[v] = u;
        }
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
    assert(r.weight == 7 && r.connected);
    assert((r.chosen == vector<int>{0, 1, 3}));

    MSTResult b = boruvka.solve();
    sort(b.chosen.begin(), b.chosen.end());
    assert(b.weight == 7 && b.connected);
    assert((b.chosen == vector<int>{0, 1, 3}));

    Kruskal forest(5);
    Boruvka boruvka_forest(5);
    for (auto [u, v, w] : vector<MSTEdge>{{0, 1, -4}, {2, 3, 6}, {3, 2, 1}}){
        forest.add_edge(u, v, w), boruvka_forest.add_edge(u, v, w);
    }
    assert(forest.solve().weight == -3 && !forest.solve().connected && forest.solve().chosen.size() == 2);
    assert(boruvka_forest.solve().weight == -3 && !boruvka_forest.solve().connected);
    assert((boruvka_forest.solve().chosen == vector<int>{0, 2}));

    Boruvka triangle(3);
    triangle.add_edge(0, 1, 7), triangle.add_edge(1, 2, 7), triangle.add_edge(2, 0, 7);
    assert(triangle.solve().weight == 14 && triangle.solve().connected && triangle.solve().chosen.size() == 2);

    assert(Kruskal(1).solve().connected && Kruskal(1).solve().weight == 0);
    assert(Kruskal(0).solve().connected && Boruvka(0).solve().connected);
    assert(!Kruskal(2).solve().connected && !Boruvka(2).solve().connected);

    const long long N = LLONG_MAX;
    DenseMSTResult d = dense_prim(vector<vector<long long>>{{N, 1, 3, 4}, {1, N, 2, N}, {3, 2, N, 5}, {4, N, 5, N}});
    assert(d.weight == 7 && d.connected);
    assert((d.parent == vector<int>{-1, 0, 1, 0}));

    vector<int> x = {5, 1, 9, 3};
    vector<vector<int>> line(4, vector<int>(4));
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++) line[i][j] = abs(x[i] - x[j]);
    }
    assert(dense_prim(line).weight == 8 && dense_prim(line).connected);

    const int M = INT_MAX;
    DenseMSTResult split = dense_prim(vector<vector<int>>{{M, M, -5}, {M, M, M}, {-5, M, M}});
    assert(split.weight == -5 && !split.connected);
    assert((split.parent == vector<int>{-1, -1, 0}));

    assert(dense_prim(vector<vector<int>>{}).connected && dense_prim(vector<vector<int>>{{0}}).weight == 0);
    return 0;
}
