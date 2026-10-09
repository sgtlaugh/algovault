/***
 *
 * Rerooting DP
 * A tree DP evaluated with every node as the root, from three user operations
 *
 * Complexity: O(n) calls to merge and O(n) calls to apply_edge, O(n) memory
 *
 * Rerooting<T, E> tree(n); tree.add_edge(u, v, w); vector<T> res = tree.solve(merge, identity, apply_edge);
 *   T: the DP value of a rooted subtree, E: the edge data (defaults to long long)
 *   apply_edge(value, from, to, w): value is the DP of the subtree hanging at from, on the far side of the
 *     edge (from, to) of weight w, the result is what that subtree contributes to to
 *   merge(a, b): combines contributions, must be associative and commutative with merge(identity, a) = a
 *   res[r] is the merge of the contributions of every neighbor of r with r as the root,
 *     identity when r has no neighbors (n = 1)
 * Node data goes inside apply_edge through from (the child's own value) or applied to res afterwards
 * Nodes are 0-based, n >= 1, solve asserts the graph is a tree, iterative so deep trees are safe
 *
 * Example: sum of distances from every node
 *   T = pair<long long, long long> (nodes below, sum of their distances), identity {0, 0}
 *   merge = {a.first + b.first, a.second + b.second}
 *   apply_edge(x, from, to, w) = {x.first + 1, x.second + (x.first + 1) * w}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename E = long long>
struct Rerooting{
    struct Edge{
        int to;
        E w;
    };

    int n, added_edges = 0;
    vector<vector<Edge>> adj;

    Rerooting(int n) : n(n), adj(n){}

    void add_edge(int u, int v, E w = E()){
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        added_edges++;
    }

    template<typename Merge, typename Apply>
    vector<T> solve(Merge merge, const T& identity, Apply apply_edge) const{
        assert(n >= 1 && added_edges == n - 1);
        vector<int> parent(n, -1), order = {0};
        parent[0] = 0;  /// the root is its own parent, so none of its neighbors is ever taken for its parent

        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            for (const Edge& e : adj[u]){
                if (e.to == parent[u]) continue;
                assert(parent[e.to] == -1);
                parent[e.to] = u;
                order.push_back(e.to);
            }
        }
        assert((int)order.size() == n);

        /// down[v]: DP of v's subtree with root 0, up[v]: DP of parent[v]'s side with v as the root
        vector<T> down(n, identity), up(n, identity), res(n, identity);
        for (int i = n - 1; i > 0; i--){
            int v = order[i];
            for (const Edge& e : adj[v]){
                if (e.to != parent[v]) down[v] = merge(down[v], apply_edge(down[e.to], e.to, v, e.w));
            }
        }

        vector<T> contribution, suffix;
        for (int v : order){
            int deg = adj[v].size();
            contribution.assign(deg, identity), suffix.assign(deg + 1, identity);
            for (int j = 0; j < deg; j++){
                const Edge& e = adj[v][j];
                contribution[j] = apply_edge(e.to == parent[v] ? up[v] : down[e.to], e.to, v, e.w);
            }
            for (int j = deg - 1; j >= 0; j--) suffix[j] = merge(contribution[j], suffix[j + 1]);
            res[v] = suffix[0];

            T prefix = identity;
            for (int j = 0; j < deg; j++){
                int c = adj[v][j].to;
                if (c != parent[v]) up[c] = merge(prefix, suffix[j + 1]);
                prefix = merge(prefix, contribution[j]);
            }
        }
        return res;
    }
};

int main(){
    /***
     *          0
     *        /   \
     *       1     2
     *      / \     \
     *     3   4     5
     *    /
     *   6
    ***/
    using Info = pair<long long, long long>;
    Rerooting<Info> tree(7);
    vector<array<int, 3>> edges = {{0, 1, 3}, {0, 2, 1}, {1, 3, 2}, {1, 4, 7}, {2, 5, 4}, {3, 6, 5}};
    for (auto [u, v, w] : edges) tree.add_edge(u, v, w);

    auto add = [](const Info& a, const Info& b){ return Info{a.first + b.first, a.second + b.second}; };
    auto lift = [](const Info& x, int, int, long long w){ return Info{x.first + 1, x.second + (x.first + 1) * w}; };
    vector<Info> sums = tree.solve(add, Info{0, 0}, lift);
    vector<long long> expected_sums = {34, 31, 37, 37, 66, 57, 62};
    for (int v = 0; v < 7; v++) assert(sums[v].first == 6 && sums[v].second == expected_sums[v]);

    /// int edge weights, one instance solved twice: weighted eccentricity, then eccentricity in edges
    Rerooting<long long, int> far(7);
    for (auto [u, v, w] : edges) far.add_edge(u, v, w);
    auto longest = [](long long a, long long b){ return max(a, b); };
    auto extend = [](long long x, int, int, int w){ return x + w; };
    auto hop = [](long long x, int, int, int){ return x + 1; };
    assert((far.solve(longest, 0LL, extend) == vector<long long>{10, 8, 11, 10, 15, 15, 15}));
    assert((far.solve(longest, 0LL, hop) == vector<long long>{3, 3, 4, 4, 4, 5, 5}));

    /***
     * Path 0 - 1 - 2 with node values 5, 1, 2: the heaviest node value among non-root nodes
    ***/
    vector<long long> value = {5, 1, 2};
    Rerooting<long long> path(3);
    path.add_edge(1, 0), path.add_edge(1, 2);
    auto heaviest = [&](long long x, int from, int, long long){ return max(x, value[from]); };
    assert((path.solve(longest, LLONG_MIN, heaviest) == vector<long long>{2, 5, 5}));

    Rerooting<long long, int> single(1);
    assert((single.solve(longest, -7LL, hop) == vector<long long>{-7}));

    return 0;
}
