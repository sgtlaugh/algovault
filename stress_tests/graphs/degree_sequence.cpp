#include "../common.h"

#define main library_main
#include "../../code_library/graphs/degree_sequence.cpp"
#undef main

/// Erdos-Gallai checked at every k without sorting tricks, O(n^2)
bool brute_graphical(vector<int> deg){
    int n = deg.size();
    long long total = 0;
    for (int d : deg){
        if (d < 0) return false;
        total += d;
    }
    if (total % 2) return false;

    sort(deg.rbegin(), deg.rend());
    for (int k = 1; k <= n; k++){
        long long lhs = 0, rhs = (long long)k * (k - 1);
        for (int i = 0; i < k; i++) lhs += deg[i];
        for (int i = k; i < n; i++) rhs += min(deg[i], k);
        if (lhs > rhs) return false;
    }

    return true;
}

int root(vector<int>& parent, int v){
    while (parent[v] != v) v = parent[v] = parent[parent[v]];

    return v;
}

void check_graph(const vector<int>& deg, const vector<pair<int, int>>& edges){
    int n = deg.size();
    vector<int> got(n, 0);
    set<pair<int, int>> seen;
    for (auto [u, v] : edges){
        assert(0 <= u && u < n && 0 <= v && v < n && u != v);
        assert(seen.insert({min(u, v), max(u, v)}).second);
        got[u]++, got[v]++;
    }
    assert(got == deg);
}

void check_tree(const vector<int>& deg, const vector<pair<int, int>>& edges){
    int n = deg.size();
    check_graph(deg, edges);
    assert((int)edges.size() == n - 1);

    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    for (auto [u, v] : edges){
        assert(root(parent, u) != root(parent, v));
        parent[root(parent, u)] = root(parent, v);
    }
}

void check_all(const vector<int>& deg, bool graphical, bool has_tree){
    assert(is_graphical(deg) == graphical);

    auto graph = graph_from_degrees(deg);
    assert(graph.has_value() == graphical);
    if (graph) check_graph(deg, *graph);

    auto tree = tree_from_degrees(deg);
    assert(tree.has_value() == has_tree);
    if (tree) check_tree(deg, *tree);
}

int main(){
    check_all({}, true, false);

    /// Every simple graph on up to 6 vertices gives the true graphical sequences, and its trees the tree sequences
    for (int n = 1; n <= 6; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }

        set<vector<int>> graphs, trees;
        int m = pairs.size();
        for (int mask = 0; mask < 1 << m; mask++){
            vector<int> deg(n, 0), parent(n);
            iota(parent.begin(), parent.end(), 0);
            int merges = 0;
            for (int e = 0; e < m; e++){
                if (!(mask >> e & 1)) continue;
                auto [u, v] = pairs[e];
                deg[u]++, deg[v]++;
                if (root(parent, u) != root(parent, v)) parent[root(parent, u)] = root(parent, v), merges++;
            }
            graphs.insert(deg);
            if (__builtin_popcount(mask) == n - 1 && merges == n - 1) trees.insert(deg);
        }

        /// Values -1 and n are out of range and must be rejected, not wrap around; n = 6 skips them for time
        int low = n < 6 ? -1 : 0, base = n < 6 ? n + 2 : n, total = 1;
        for (int i = 0; i < n; i++) total *= base;
        for (int code = 0; code < total; code++){
            vector<int> deg(n);
            for (int i = 0, c = code; i < n; i++, c /= base) deg[i] = c % base + low;
            check_all(deg, graphs.count(deg), trees.count(deg));
        }
    }

    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 120);
        double p = stress::rand_int(0, 100) / 100.0;
        vector<int> deg(n, 0);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (stress::rand_int(0, 999) < p * 1000) deg[u]++, deg[v]++;
            }
        }

        /// Shifting one unit between two vertices often crosses the graphical boundary
        if (it % 2){
            int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 1);
            deg[a]++, deg[b]--;
        }
        if (it % 5 == 0) deg[stress::rand_int(0, n - 1)] += stress::rand_int(-2, 2);

        bool graphical = brute_graphical(deg);
        assert(is_graphical(deg) == graphical);
        auto graph = graph_from_degrees(deg);
        assert(graph.has_value() == graphical);
        if (graph) check_graph(deg, *graph);
    }

    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(2, 2000);
        vector<int> deg(n, 0);
        for (int v = 1; v < n; v++){
            int parent = it % 2 ? stress::rand_int(max(0, v - 3), v - 1) : stress::rand_int(0, v - 1);
            deg[v]++, deg[parent]++;
        }
        shuffle(deg.begin(), deg.end(), stress::rng());
        check_tree(deg, *tree_from_degrees(deg));
        assert(is_graphical(deg));

        int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 1);
        if (a != b && deg[a] > 1){
            deg[a]--, deg[b]++;
            check_tree(deg, *tree_from_degrees(deg));
            deg[b]++;
            assert(!tree_from_degrees(deg));
        }
    }

    /// The complete graph on 700 vertices, about 2.4e5 edges, and a star
    int n = 700;
    vector<int> full(n, n - 1);
    check_graph(full, *graph_from_degrees(full));
    full[0]++;
    assert(!is_graphical(full));

    vector<int> star(n, 1);
    star[0] = n - 1;
    check_graph(star, *graph_from_degrees(star));
    check_tree(star, *tree_from_degrees(star));

    return 0;
}
