#include "../common.h"

#define main library_main
#include "../../code_library/graphs/minimum_spanning_tree.cpp"
#undef main

/// Union by size with a two-pass find, written apart from MSTForest so the checks do not trust the code under test
struct ReferenceForest{
    vector<int> parent, size;

    ReferenceForest(int n) : parent(n), size(n, 1) {
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int x){
        int root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != root){
            int next = parent[x];
            parent[x] = root, x = next;
        }
        return root;
    }

    bool join(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;

        if (size[a] < size[b]) swap(a, b);
        parent[b] = a, size[a] += size[b];
        return true;
    }
};

/// Prim's algorithm on an adjacency matrix, run from every unvisited vertex so it covers each component
pair<long long, int> prim(int n, const vector<MSTEdge>& edges){
    const long long INF = LLONG_MAX;
    vector<vector<long long>> w(n, vector<long long>(n, INF));
    for (auto [u, v, c] : edges){
        if (u != v) w[u][v] = w[v][u] = min(w[u][v], c);
    }

    vector<char> used(n, 0);
    vector<long long> best(n, INF);
    long long total = 0;
    int components = 0;

    for (int it = 0; it < n; it++){
        int u = -1;
        for (int v = 0; v < n; v++){
            if (!used[v] && (u == -1 || best[v] < best[u])) u = v;
        }
        if (best[u] == INF) components++;
        else total += best[u];
        used[u] = 1;
        for (int v = 0; v < n; v++){
            if (!used[v] && w[u][v] < best[v]) best[v] = w[u][v];
        }
    }
    return {total, components};
}

/// Every edge subset: a minimum spanning forest is a lightest acyclic subset among those with the most edges
pair<long long, int> brute(int n, const vector<MSTEdge>& edges){
    int m = edges.size(), most = -1;
    long long lightest = 0;

    for (int mask = 0; mask < (1 << m); mask++){
        ReferenceForest dsu(n);
        bool acyclic = true;
        long long sum = 0;
        for (int i = 0; i < m && acyclic; i++){
            if (!(mask >> i & 1)) continue;
            acyclic = dsu.join(edges[i].u, edges[i].v);
            sum += edges[i].w;
        }

        int size = __builtin_popcount(mask);
        if (!acyclic) continue;
        if (size > most || (size == most && sum < lightest)) most = size, lightest = sum;
    }
    return {lightest, n - most};
}

void check_edges(int n, const vector<MSTEdge>& edges, const MSTResult& r, pair<long long, int> expected){
    auto [total, components] = expected;
    assert(r.weight == total);
    assert((int)r.chosen.size() == n - components);
    assert(r.connected == (components <= 1));

    /// The chosen edges must be distinct and form a forest: union-find never sees a cycle
    ReferenceForest dsu(n);
    long long sum = 0;
    for (int i : r.chosen){
        assert(dsu.join(edges[i].u, edges[i].v));
        sum += edges[i].w;
    }
    assert(sum == r.weight);
}

template<typename T>
void check_dense(const vector<vector<T>>& w, pair<long long, int> expected){
    int n = w.size();
    DenseMSTResult r = dense_prim(w);
    auto [total, components] = expected;
    assert(r.weight == total);
    assert(r.connected == (components <= 1));

    /// parent edges must exist, sum to the weight and form a forest with one root per component
    ReferenceForest dsu(n);
    long long sum = 0;
    int roots = 0;
    for (int v = 0; v < n; v++){
        if (r.parent[v] == -1){
            roots++;
            continue;
        }
        assert(r.parent[v] != v && w[v][r.parent[v]] != numeric_limits<T>::max());
        assert(dsu.join(v, r.parent[v]));
        sum += w[v][r.parent[v]];
    }
    assert(sum == r.weight && roots == components);
}

void check(int n, const vector<MSTEdge>& edges, pair<long long, int> expected, bool dense){
    Kruskal kruskal(n);
    Boruvka boruvka(n);
    for (auto [u, v, w] : edges) kruskal.add_edge(u, v, w), boruvka.add_edge(u, v, w);
    check_edges(n, edges, kruskal.solve(), expected);
    check_edges(n, edges, boruvka.solve(), expected);
    if (!dense) return;

    vector<vector<long long>> w(n, vector<long long>(n, LLONG_MAX));
    for (auto [u, v, c] : edges){
        if (u != v) w[u][v] = w[v][u] = min(w[u][v], c);
    }
    check_dense(w, expected);
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 7), m = n ? stress::rand_int(0, 11) : 0;
        long long range = it % 3 == 0 ? 2 : 1000000000000LL;
        vector<MSTEdge> edges;
        for (int i = 0; i < m; i++){
            edges.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(-range, range)});
        }

        auto expected = brute(n, edges);
        assert(expected == prim(n, edges));
        check(n, edges, expected, true);
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, 40), m = n ? stress::rand_int(0, 3 * n) : 0;
        long long range = it % 3 == 0 ? 3 : 1000000000000LL;
        vector<MSTEdge> edges;
        for (int i = 0; i < m; i++){
            edges.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(-range, range)});
        }
        check(n, edges, prim(n, edges), true);
    }

    /// A large instance with a known answer: tree edge (parent[i] < i, i) weighs i * K - C, so the heaviest tree edge
    /// on any u-v path ends at max(u, v) and an extra edge weighing at least max(u, v) * K - C can never improve the tree
    /// The last 3 vertices stay isolated, and vertices are relabeled so indices carry no order
    const long long K = 1000, C = 100000000;
    int n = 200000, isolated = 3, core = n - isolated;
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<MSTEdge> edges;
    long long total = 0;
    for (int i = 1; i < core; i++){
        int p = i % 3 ? i - 1 : stress::rand_int(0, i - 1);
        edges.push_back({label[p], label[i], i * K - C});
        total += i * K - C;
    }
    for (int i = 0; i < n; i++){
        int u = stress::rand_int(0, core - 1), v = stress::rand_int(0, core - 1);
        edges.push_back({label[u], label[v], max(u, v) * K - C + stress::rand_int(0, i % 2 ? 0 : 50 * K)});
    }

    shuffle(edges.begin(), edges.end(), stress::rng());
    check(n, edges, {total, 1 + isolated}, false);

    /// A complete graph with int weights in a narrow range, so ties are everywhere: Boruvka and dense_prim against Kruskal
    int dn = 1500;
    vector<vector<int>> w(dn, vector<int>(dn, 0));
    Kruskal complete(dn);
    Boruvka complete_boruvka(dn);
    for (int u = 0; u < dn; u++){
        for (int v = u + 1; v < dn; v++){
            w[u][v] = w[v][u] = stress::rand_int(-20, 20);
            complete.add_edge(u, v, w[u][v]), complete_boruvka.add_edge(u, v, w[u][v]);
        }
    }
    long long expected = complete.solve().weight;
    assert(complete_boruvka.solve().weight == expected);
    check_dense(w, {expected, 1});
    return 0;
}
