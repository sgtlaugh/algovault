/***
 *
 * Rerooting DP
 * A tree DP evaluated with every node as the root, from three user operations
 *
 * Complexity: O(n), 2(n - 1) calls to apply_edge and O(n) calls to combine
 *
 * Rerooting<T, E = long long> tree(n): nodes 0..n-1, T is the DP value of a rooted subtree and E the edge data
 * tree.add_edge(u, v, w = E()): adds the edge (u, v) with data w, 0 <= u, v < n
 * tree.solve(combine, identity, apply_edge): res with res[r] the contributions of every neighbor of r combined
 *
 * Notes:
 *   - apply_edge(value, from, to, w): what the subtree at from, whose DP is value, contributes to to across edge w
 *   - combine(a, b): two contributions as one, associative and commutative with combine(identity, a) = a
 *   - solve asserts n >= 1 and a tree, res[0] is identity when n = 1
 *   - Node data goes inside apply_edge through from, the child's own node, or is applied to res afterwards
 *   - E must be default constructible
 *   - Iterative, so deep trees are safe
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<typename T, typename E = long long>
struct Rerooting{
    int n;
    vector<array<int, 2>> ends;
    vector<E> data;

    Rerooting(int n) : n(n) {}

    void add_edge(int u, int v, E w = E()){
        assert(0 <= u && u < n && 0 <= v && v < n);
        ends.push_back({u, v});
        data.push_back(w);
    }

    template<typename Combine, typename Apply>
    vector<T> solve(Combine combine, const T& identity, Apply apply_edge) const{
        assert(n >= 1 && (int)ends.size() == n - 1);

        /// The edges at v are adj[start[v]..start[v + 1] - 1] in insertion order, as indices into ends and data
        vector<int> start(n + 1, 0), adj(2 * (n - 1));
        for (auto [u, v] : ends) start[u + 1]++, start[v + 1]++;
        for (int v = 0; v < n; v++) start[v + 1] += start[v];
        vector<int> pos(start.begin(), start.end() - 1);
        for (int e = 0; e < n - 1; e++) adj[pos[ends[e][0]]++] = e, adj[pos[ends[e][1]]++] = e;

        /// Indexed by BFS position from 0, so the children of i are first[i]..first[i + 1] - 1 and passes scan in order
        vector<int> order = {0}, parent(n, -1), edge = {-1}, first(n + 1, n);
        vector<char> seen(n, false);
        seen[0] = true;
        for (int i = 0; i < (int)order.size(); i++){
            int u = order[i];
            first[i] = order.size();
            for (int j = start[u]; j < start[u + 1]; j++){
                int e = adj[j], to = ends[e][0] ^ ends[e][1] ^ u;
                if (e == edge[i]) continue;
                assert(!seen[to]);
                seen[to] = true, parent[order.size()] = i;
                order.push_back(to), edge.push_back(e);
            }
        }
        assert((int)order.size() == n);

        /// down[i] combines i's children, then what i's subtree gives its parent, up[i] what the parent side gives i
        vector<T> down(n, identity), up(n, identity), suffix, res(n, identity);
        for (int i = n - 1; i > 0; i--){
            down[i] = apply_edge(down[i], order[i], order[parent[i]], data[edge[i]]);
            down[parent[i]] = combine(down[parent[i]], down[i]);
        }

        for (int i = 0; i < n; i++){
            int l = first[i], k = first[i + 1] - l;
            suffix.assign(k + 1, identity);
            for (int j = k - 1; j >= 0; j--) suffix[j] = combine(down[l + j], suffix[j + 1]);
            res[order[i]] = combine(up[i], suffix[0]);

            T prefix = up[i];
            for (int j = 0; j < k; j++){
                up[l + j] = apply_edge(combine(prefix, suffix[j + 1]), order[i], order[l + j], data[edge[l + j]]);
                prefix = combine(prefix, down[l + j]);
            }
        }

        return res;
    }
};

int main(){
    /***
     *             0
     *        3 /     \ 1
     *         1       2
     *    2 /    \ 7    \ 4
     *     3      4      5
     *  5 /
     *   6
    ***/
    vector<array<int, 3>> edges = {{0, 1, 3}, {0, 2, 1}, {1, 3, 2}, {1, 4, 7}, {2, 5, 4}, {3, 6, 5}};

    /// Sum of distances from every node, a value is {nodes below, sum of their distances}
    using Info = pair<long long, long long>;
    Rerooting<Info> tree(7);
    for (auto [u, v, w] : edges) tree.add_edge(u, v, w);
    auto add = [](const Info& a, const Info& b){ return Info{a.first + b.first, a.second + b.second}; };
    auto lift = [](const Info& x, int, int, long long w){ return Info{x.first + 1, x.second + (x.first + 1) * w}; };
    vector<Info> sums = tree.solve(add, Info{0, 0}, lift);
    assert(sums[0].second == 34);  /// 3 + 1 + 5 + 10 + 5 + 10
    assert(sums[6].second == 62);  /// 5 + 7 + 10 + 11 + 14 + 15

    /// Eccentricity, the farthest node from each root, with int edge weights
    Rerooting<long long, int> far(7);
    for (auto [u, v, w] : edges) far.add_edge(u, v, w);
    auto longest = [](long long a, long long b){ return max(a, b); };
    auto extend = [](long long x, int, int, int w){ return x + w; };
    assert((far.solve(longest, 0LL, extend) == vector<long long>{10, 8, 11, 10, 15, 15, 15}));

    return 0;
}
