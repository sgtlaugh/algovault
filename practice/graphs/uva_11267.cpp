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
 *     which node of a tree becomes its root is unspecified
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

    MSTForest(int n) : parent(n), rank(n, 0){
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
    int n;
    while (scanf("%d", &n) == 1 && n){
        int e;
        scanf("%d", &e);

        Kruskal mst(n);
        vector<vector<int>> adj(n);
        for (int i = 0; i < e; i++){
            int a, b, c;
            scanf("%d %d %d", &a, &b, &c);
            mst.add_edge(a - 1, b - 1, c);
            adj[a - 1].push_back(b - 1), adj[b - 1].push_back(a - 1);
        }

        /// Every pair is a buyer and a coder, so the network must be two-colourable
        vector<int> side(n, -1), queue = {0};
        side[0] = 0;
        bool bipartite = true;
        for (int i = 0; i < (int)queue.size() && bipartite; i++){
            int u = queue[i];
            for (int v: adj[u]){
                if (side[v] == -1) side[v] = side[u] ^ 1, queue.push_back(v);
                else if (side[v] == side[u]) bipartite = false;
            }
        }
        if (!bipartite){
            puts("Invalid data, Idiot!");
            continue;
        }

        /// The tree keeps everyone connected, every other profitable pairing only lowers the loss
        MSTResult r = mst.solve();
        vector<bool> in_tree(e, false);
        for (int i: r.chosen) in_tree[i] = true;

        long long loss = r.weight;
        for (int i = 0; i < e; i++){
            if (!in_tree[i] && mst.edges[i].w < 0) loss += mst.edges[i].w;
        }
        printf("%lld\n", loss);
    }
    return 0;
}
