/***
 *
 * Centroid Decomposition
 * Recursively cuts a weighted tree (or forest) at centroids, the path u-v always passes through lca(u, v) in the centroid tree
 *
 * Complexity: O(n log n) to build, plus the visitor's own work on O(n log n) delivered (node, distance) pairs
 *
 * CentroidDecomposition cd(n); cd.add_edge(u, v, w); cd.build(visit);  (w defaults to 1, nodes 0-based)
 *   cd.parent[c]: parent of c in the centroid tree, -1 for the root centroid of each tree of the forest
 *   cd.depth[c]: depth of c in the centroid tree, at most floor(log2(n))
 *   visit(c, branches) is called once per centroid, parents before children, with c already cut out:
 *     branches[i] lists (node, distance to c) for every node of the i-th piece the cut of c leaves,
 *     the piece hanging off one neighbor of c, which becomes the component of one child centroid
 *     its type is CentroidDecomposition::Branches, a vector<vector<pair<int, long long>>>
 *   build() without a visitor only fills parent and depth
 * Distances are long long sums of edge weights, keep 2 * n * max|w| below 9.2e18 so d1 + d2 cannot overflow
 * Iterative, so deep trees cannot overflow the stack
 *
 * Example, count_paths_at_most below: unordered pairs u != v with dist(u, v) <= k in O(n log^2 n)
 *   at each centroid, count pairs over {c} and all branches with d1 + d2 <= k (sort + two pointers),
 *   minus the pairs inside a single branch, whose path does not pass through c
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct CentroidDecomposition{
    using Branches = vector<vector<pair<int, long long>>>;

    int n;
    vector<vector<pair<int, long long>>> adj;
    vector<int> parent, depth;

    CentroidDecomposition(int n) : n(n), adj(n), parent(n, -1), depth(n, 0), removed(n, 0), from(n, -1), size(n, 0) {}

    void add_edge(int u, int v, long long w = 1){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    /// Skips collecting branches, which is about a third of the build time
    void build(){
        decompose([](int){});
    }

    template <typename Visit>
    void build(Visit visit){
        decompose([&](int c){ visit(c, collect_branches(c)); });
    }

private:
    vector<char> removed;
    vector<int> from, size, order;

    Branches collect_branches(int c){
        Branches branches;
        for (auto [v, w] : adj[c]){
            if (removed[v]) continue;
            vector<pair<int, long long>> branch = {{v, w}};
            from[v] = c;
            for (int i = 0; i < (int)branch.size(); i++){
                auto [x, d] = branch[i];  /// a copy: pushing below may reallocate the branch
                for (auto [y, wy] : adj[x]){
                    if (y != from[x] && !removed[y]) from[y] = x, branch.push_back({y, d + wy});
                }
            }
            branches.push_back(move(branch));
        }
        return branches;
    }

    template <typename OnCentroid>
    void decompose(OnCentroid on_centroid){
        removed.assign(n, 0);

        /// (any node of a component, centroid of the component it was cut from)
        vector<pair<int, int>> stack;
        for (int s = 0; s < n; s++){
            if (!removed[s]) stack.push_back({s, -1});
            while (!stack.empty()){
                auto [r, up] = stack.back();
                stack.pop_back();
                int c = find_centroid(r);
                removed[c] = 1, parent[c] = up, depth[c] = up == -1 ? 0 : depth[up] + 1;

                on_centroid(c);
                for (auto [v, w] : adj[c]){
                    if (!removed[v]) stack.push_back({v, c});
                }
            }
        }
    }

    /// Walks from r toward the child holding more than half the component, the part behind stays below half
    int find_centroid(int r){
        order.assign(1, r);
        from[r] = -1;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (auto [v, w] : adj[u]){
                if (v != from[u] && !removed[v]) from[v] = u, order.push_back(v);
            }
        }

        int total = order.size();
        for (int i = total - 1; i >= 0; i--){
            int u = order[i];
            size[u] = 1;
            for (auto [v, w] : adj[u]){
                if (v != from[u] && !removed[v]) size[u] += size[v];
            }
        }

        for (int u = r;;){
            int next = -1;
            for (auto [v, w] : adj[u]){
                if (v != from[u] && !removed[v] && size[v] > total / 2) next = v;
            }
            if (next == -1) return u;
            u = next;
        }
    }
};

long long count_pairs_at_most(vector<long long> d, long long k){
    sort(d.begin(), d.end());

    long long pairs = 0;
    int j = (int)d.size() - 1;
    for (int i = 0; i < j; i++){
        while (i < j && d[i] + d[j] > k) j--;
        pairs += j - i;
    }
    return pairs;
}

long long count_paths_at_most(CentroidDecomposition& cd, long long k){
    long long pairs = 0;
    cd.build([&](int, const CentroidDecomposition::Branches& branches){
        vector<long long> all = {0};
        for (auto& branch : branches){
            vector<long long> d;
            for (auto [v, x] : branch) d.push_back(x);
            pairs -= count_pairs_at_most(d, k);
            all.insert(all.end(), d.begin(), d.end());
        }
        pairs += count_pairs_at_most(all, k);
    });
    return pairs;
}

int main(){
    /// Path 0-1-2-3-4: distance 1 for 4 pairs, 2 for 3, 3 for 2, 4 for 1
    CentroidDecomposition path(5);
    for (int i = 1; i < 5; i++) path.add_edge(i - 1, i);
    assert(count_paths_at_most(path, 0) == 0);
    assert(count_paths_at_most(path, 1) == 4);
    assert(count_paths_at_most(path, 2) == 7);
    assert(count_paths_at_most(path, 4) == 10);
    assert((path.parent == vector<int>{1, 2, -1, 2, 3}));
    assert((path.depth == vector<int>{2, 1, 0, 1, 2}));

    /// Star at 0 with weights 1, 2, 3: distances 1, 2, 3 from the center, 3, 4, 5 between leaves
    CentroidDecomposition star(4);
    for (int i = 1; i < 4; i++) star.add_edge(0, i, i);
    assert(count_paths_at_most(star, 3) == 4);
    assert(count_paths_at_most(star, 4) == 5);
    assert(count_paths_at_most(star, 5) == 6);

    vector<int> visited;
    star.build([&](int c, const CentroidDecomposition::Branches& branches){
        visited.push_back(c);
        if (c != 0) assert(branches.empty());
        else assert((branches == CentroidDecomposition::Branches{{{1, 1}}, {{2, 2}}, {{3, 3}}}));
    });
    assert(visited.size() == 4 && visited[0] == 0);
    assert((star.parent == vector<int>{-1, 0, 0, 0}));

    CentroidDecomposition forest(4);
    forest.add_edge(0, 1), forest.add_edge(2, 3);
    assert(count_paths_at_most(forest, 1) == 2);
    assert(count_paths_at_most(forest, 100) == 2);
    assert(count(forest.parent.begin(), forest.parent.end(), -1) == 2);

    CentroidDecomposition single(1);
    single.build();
    assert(single.parent[0] == -1 && single.depth[0] == 0);
    assert(count_paths_at_most(single, 5) == 0);

    CentroidDecomposition empty(0);
    empty.build();
    assert(count_paths_at_most(empty, 3) == 0);

    return 0;
}
