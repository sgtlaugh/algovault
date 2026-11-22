/***
 *
 * Degree Sequence
 * Whether a list of degrees belongs to a simple graph, and a simple graph or tree that has them
 *
 * Complexity:
 *   - is_graphical: O(n log n)
 *   - graph_from_degrees: O((n + m) log n) for m = sum(deg) / 2 edges
 *   - tree_from_degrees: O(n)
 *
 * graph_from_degrees(deg): the edges of a simple graph where vertex i has degree deg[i], nullopt if none exists
 * is_graphical(deg): true if some simple graph has vertex i of degree deg[i], for vertices 0 to n - 1
 * tree_from_degrees(deg): the n - 1 edges of a tree where vertex i has degree deg[i], nullopt if none exists
 *
 * Notes:
 *   - Simple means no self loops and no multi-edges; any int values are allowed, a negative one is never graphical
 *   - A tree exists iff n = 1 with deg = {0}, or n >= 2 with every degree >= 1 and the sum 2 (n - 1)
 *   - An empty deg is graphical with no edges, but has no tree
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Erdos-Gallai: with deg sorted descending, the k largest sum to at most k (k - 1) + sum of min(deg[i], k) over i >= k
bool is_graphical(vector<int> deg){
    int n = deg.size();
    long long total = 0;
    for (int d : deg){
        if (d < 0) return false;
        total += d;
    }
    if (total % 2) return false;

    sort(deg.rbegin(), deg.rend());
    vector<long long> suffix(n + 1, 0);
    for (int i = n - 1; i >= 0; i--) suffix[i] = suffix[i + 1] + deg[i];

    long long prefix = 0;
    int at_least_k = n;
    for (int k = 1; k <= n; k++){
        prefix += deg[k - 1];
        while (at_least_k > 0 && deg[at_least_k - 1] < k) at_least_k--;

        int split = max(k, at_least_k);
        if (prefix > (long long)k * (k - 1) + (long long)k * (split - k) + suffix[split]) return false;
    }

    return true;
}

optional<vector<pair<int, int>>> graph_from_degrees(const vector<int>& deg){
    if (!is_graphical(deg)) return nullopt;

    /// Havel-Hakimi: joining the largest remaining degree to the next largest ones keeps the rest graphical
    set<pair<int, int>> open;
    for (int v = 0; v < (int)deg.size(); v++){
        if (deg[v] > 0) open.insert({deg[v], v});
    }

    vector<pair<int, int>> edges, touched;
    while (!open.empty()){
        auto [d, u] = *open.rbegin();
        open.erase(prev(open.end()));

        touched.clear();
        for (int i = 0; i < d; i++){
            auto [dv, v] = *open.rbegin();
            open.erase(prev(open.end()));
            edges.push_back({u, v});
            touched.push_back({dv - 1, v});
        }
        for (auto p : touched){
            if (p.first > 0) open.insert(p);
        }
    }

    return edges;
}

optional<vector<pair<int, int>>> tree_from_degrees(const vector<int>& deg){
    int n = deg.size();
    if (n == 1 && deg[0] == 0) return vector<pair<int, int>>{};
    if (n < 2) return nullopt;

    long long total = 0;
    for (int d : deg){
        if (d < 1) return nullopt;
        total += d;
    }
    if (total != 2LL * (n - 1)) return nullopt;

    vector<int> leaves, inner;
    for (int v = 0; v < n; v++) (deg[v] == 1 ? leaves : inner).push_back(v);
    if (inner.empty()) return vector<pair<int, int>>{{0, 1}};

    /// A caterpillar: the inner vertices form a path and their free slots add up to exactly the number of leaves
    vector<pair<int, int>> edges;
    int next_leaf = 0, k = inner.size();
    for (int i = 0; i < k; i++){
        int u = inner[i], slots = deg[u] - (i > 0) - (i + 1 < k);
        if (i > 0) edges.push_back({inner[i - 1], u});
        for (int j = 0; j < slots; j++) edges.push_back({u, leaves[next_leaf++]});
    }

    return edges;
}

int main(){
    /// 3, 3, 3, 3 is the complete graph on 4 vertices; in 3, 3, 1, 1 both 3s touch every vertex, so no 1 stays 1
    assert(is_graphical({3, 3, 3, 3}));
    assert(!is_graphical({3, 3, 1, 1}));
    assert(!is_graphical({2, 1, 0}));
    assert(!is_graphical({4, 1, 1, 1}));

    /// Havel-Hakimi on 2, 2, 1, 1: vertex 1 joins 0 and 3, leaving 1, 0, 1, 0, then 2 joins 0
    vector<pair<int, int>> path = *graph_from_degrees({2, 2, 1, 1});
    assert((path == vector<pair<int, int>>{{1, 0}, {1, 3}, {2, 0}}));
    assert(!graph_from_degrees({1, 1, 1}));

    /// Leaves 0, 2, 3, 5 and inner path 1 - 4: vertex 1 takes leaves 0 and 2, vertex 4 takes 3 and 5
    vector<pair<int, int>> tree = *tree_from_degrees({1, 3, 1, 1, 3, 1});
    assert((tree == vector<pair<int, int>>{{1, 0}, {1, 2}, {1, 4}, {4, 3}, {4, 5}}));
    assert(!tree_from_degrees({2, 2, 2}));
    assert(tree_from_degrees({0})->empty());

    return 0;
}
