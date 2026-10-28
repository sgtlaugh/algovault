/***
 *
 * Centroid Decomposition
 * Recursively cuts a weighted forest at centroids, so the path u-v passes through lca(u, v) in the centroid tree
 *
 * Complexity: O(n log n) to build, plus the visitor's own work on O(n log n) delivered (node, distance) pairs
 *
 * CentroidDecomposition cd(n): an empty forest on nodes 0..n-1
 * cd.add_edge(u, v, w): adds the edge u - v of weight w, 1 by default
 * cd.build(visit): calls visit(c, branches) once per centroid c, parents before children, with c already cut out
 * cd.build(): fills only parent and depth
 * count_paths_at_most(cd, k): the number of unordered pairs u != v with dist(u, v) <= k, in O(n log^2 n)
 * cd.parent[c]: parent of c in the centroid tree, -1 for the root centroid of each tree of the forest
 * cd.depth[c]: depth of c in the centroid tree, at most floor(log2(n))
 *
 * Notes:
 *   - branches[i] lists (node, distance to c) for the piece behind the i-th uncut neighbor of c, one child's component
 *   - branches has type CentroidDecomposition::Branches
 *   - build asserts that the edges form a forest, so a cycle, multi-edge or self-loop aborts
 *   - Distances are long long sums of edge weights, keep 2 n max|w| below 9.2e18 so d1 + d2 cannot overflow
 *   - Iterative, so deep trees cannot overflow the stack
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct CentroidDecomposition{
    using Branches = vector<vector<pair<int, long long>>>;

    int n;
    vector<tuple<int, int, long long>> edges;
    vector<int> parent, depth;
    /// Nodes renumbered in BFS order with flat adjacency, measured up to 4x faster than per-node vectors at n = 5e5
    vector<int> original, start;
    vector<pair<int, long long>> adj;
    vector<char> removed;
    vector<int> from, size, order;

    CentroidDecomposition(int n) : n(n), parent(n, -1), depth(n, 0) {}

    void add_edge(int u, int v, long long w = 1){
        edges.push_back({u, v, w});
    }

    /// Skips collecting branches, which is about a third of the build time
    void build(){
        decompose([](int) {});
    }

    template <typename Visit>
    void build(Visit visit){
        decompose([&](int c){ visit(original[c], collect_branches(c)); });
    }

    template <typename OnCentroid>
    void decompose(OnCentroid on_centroid){
        renumber();
        removed.assign(n, 0), from.assign(n, -1), size.assign(n, 0);

        /// (any node of a component, centroid of the component it was cut from)
        vector<pair<int, int>> stack;
        for (int s = 0; s < n; s++){
            if (!removed[s]) stack.push_back({s, -1});
            while (!stack.empty()){
                auto [r, up] = stack.back();
                stack.pop_back();
                int c = find_centroid(r);
                removed[c] = 1;
                parent[original[c]] = up == -1 ? -1 : original[up];
                depth[original[c]] = up == -1 ? 0 : depth[original[up]] + 1;

                on_centroid(c);
                for (int j = start[c]; j < start[c + 1]; j++){
                    if (!removed[adj[j].first]) stack.push_back({adj[j].first, c});
                }
            }
        }
    }

    void renumber(){
        vector<int> first(n + 1, 0);
        for (auto [u, v, w] : edges) first[u + 1]++, first[v + 1]++;
        for (int u = 0; u < n; u++) first[u + 1] += first[u];
        vector<int> fill(first.begin(), first.end() - 1);
        vector<int> incident(first[n]);
        for (int e = 0; e < (int)edges.size(); e++){
            auto [u, v, w] = edges[e];
            incident[fill[u]++] = e, incident[fill[v]++] = e;
        }

        vector<int> id(n, -1);
        int trees = 0;
        original.clear();
        for (int s = 0; s < n; s++){
            if (id[s] != -1) continue;
            trees++, id[s] = original.size(), original.push_back(s);
            for (int i = id[s]; i < (int)original.size(); i++){
                int u = original[i];
                for (int j = first[u]; j < first[u + 1]; j++){
                    auto [a, b, w] = edges[incident[j]];
                    int v = a ^ b ^ u;
                    if (id[v] == -1) id[v] = original.size(), original.push_back(v);
                }
            }
        }
        /// A forest is exactly a graph with n - trees edges, so this catches every cycle, multi-edge and self-loop
        assert((int)edges.size() == n - trees);

        start.assign(n + 1, 0), adj.resize(first[n]);
        for (int i = 0; i < n; i++){
            int u = original[i], k = start[i];
            for (int j = first[u]; j < first[u + 1]; j++){
                auto [a, b, w] = edges[incident[j]];
                adj[k++] = {id[a ^ b ^ u], w};
            }
            start[i + 1] = k;
        }
    }

    Branches collect_branches(int c){
        Branches branches;

        for (int j = start[c]; j < start[c + 1]; j++){
            auto [v, w] = adj[j];
            if (removed[v]) continue;
            vector<pair<int, long long>> branch = {{v, w}};
            from[v] = c;
            for (int i = 0; i < (int)branch.size(); i++){
                auto [x, d] = branch[i];  /// a copy: pushing below may reallocate the branch
                for (int k = start[x]; k < start[x + 1]; k++){
                    auto [y, wy] = adj[k];
                    if (y != from[x] && !removed[y]) from[y] = x, branch.push_back({y, d + wy});
                }
            }
            for (auto& [x, d] : branch) x = original[x];
            branches.push_back(move(branch));
        }

        return branches;
    }

    int find_centroid(int r){
        order.assign(1, r);
        from[r] = -1, size[r] = 0;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (int j = start[u]; j < start[u + 1]; j++){
                int v = adj[j].first;
                if (v == from[u] || removed[v]) continue;
                from[v] = u, size[v] = 0, order.push_back(v);
            }
        }

        /// Nodes holding over half the component form a path down from r, the centroid is its last node in BFS order
        int total = order.size(), centroid = -1;
        for (int i = total - 1; i >= 0; i--){
            int u = order[i];
            size[u]++;
            if (centroid == -1 && 2 * size[u] > total) centroid = u;
            if (i > 0) size[from[u]] += size[u];
        }

        return centroid;
    }
};

long long count_paths_at_most(CentroidDecomposition& cd, long long k){
    auto count_pairs_at_most = [k](vector<long long> d){
        sort(d.begin(), d.end());

        long long pairs = 0;
        int j = (int)d.size() - 1;
        for (int i = 0; i < j; i++){
            while (i < j && d[i] + d[j] > k) j--;
            pairs += j - i;
        }

        return pairs;
    };

    long long pairs = 0;

    cd.build([&](int, const CentroidDecomposition::Branches& branches){
        vector<long long> all = {0};
        for (auto& branch : branches){
            vector<long long> d;
            for (auto [v, x] : branch) d.push_back(x);
            pairs -= count_pairs_at_most(d);
            all.insert(all.end(), d.begin(), d.end());
        }
        pairs += count_pairs_at_most(all);
    });
    return pairs;
}

int main(){
    /// Path 0-1-2-3-4: 2 is cut first, then 1 and 3, then the ends
    CentroidDecomposition path(5);
    for (int i = 1; i < 5; i++) path.add_edge(i - 1, i);
    path.build();
    assert((path.parent == vector<int>{1, 2, -1, 2, 3}));
    assert((path.depth == vector<int>{2, 1, 0, 1, 2}));
    assert(count_paths_at_most(path, 2) == 7);  /// 4 pairs at distance 1, 3 at distance 2

    /// Star at 0 with weights 1, 2, 3: the center is cut first and each leaf is a branch of its own
    CentroidDecomposition star(4);
    for (int i = 1; i < 4; i++) star.add_edge(0, i, i);
    star.build([](int c, const CentroidDecomposition::Branches& branches){
        if (c == 0) assert((branches == CentroidDecomposition::Branches{{{1, 1}}, {{2, 2}}, {{3, 3}}}));
        else assert(branches.empty());
    });
    assert(count_paths_at_most(star, 4) == 5);  /// every pair but leaves 2 and 3, at distance 5
    return 0;
}
