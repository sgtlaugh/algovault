/***
 *
 * Mo's Algorithm On Tree Paths
 * Answers offline path queries on a tree by running Mo's algorithm over an Euler tour that lists each node twice
 *
 * Complexity: O(n log n) to build, O(n * sqrt(q) + q) calls to add / remove per run, plus O(q log q) sorting
 *
 * TreePathMo tree(n); tree.add_edge(u, v); tree.build(root);
 * tree.run(paths, add, remove, get, edges) returns the answers, res[qi] = get(qi):
 *   paths[i] = {u, v}, add(node) / remove(node) receive tree nodes, get(qi) returns the answer for the current path
 *   edges = false: the window holds every node on the path u..v
 *   edges = true: the window holds every node on the path except the lca, store the value of edge (parent[v], v) on v
 *   The window starts and ends empty, so the same callbacks and state can serve several runs
 *
 * Embeds a checked copy of mo.cpp's mo to stay standalone
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN COPY mo from code_library/data_structures/mo.cpp
template <typename Add, typename Remove, typename Get>
auto mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Get get){
    int q = queries.size();
    vector<decay_t<decltype(get(0))>> res(q);
    if (q == 0) return res;
    int block = max(1, (int)(n / sqrt((double)q)));

    /// Odd-even block order: Hilbert curve order measured about 5% more add / remove calls on random n = q = 5e5
    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? queries[a].second > queries[b].second : queries[a].second < queries[b].second;
    });

    /// Growing before shrinking keeps every removed position inside the window, empty targets included
    int cur_l = 0, cur_r = -1;
    for (int qi : order){
        auto [l, r] = queries[qi];
        assert(0 <= l && l <= r + 1 && r < n);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        res[qi] = get(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);  /// leave the window empty so the caller's state can be reused
    return res;
}
/// END COPY mo

struct TreePathMo{
    int n, lg;
    vector<vector<int>> adj, up;
    vector<int> depth, first, last, tour;

    TreePathMo(int n) : n(n), adj(n), depth(n), first(n), last(n){
        lg = 1;
        while ((1 << lg) < n) lg++;
    }

    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    /// Euler tour with every node at its entry and exit, iterative so deep trees are safe
    void build(int root = 0){
        up.assign(lg + 1, vector<int>(n, root));
        tour.clear();
        vector<int> it(n, 0), stack = {root};
        depth[root] = 0, first[root] = 0, tour.push_back(root);

        while (!stack.empty()){
            int u = stack.back();
            if (it[u] < (int)adj[u].size()){
                int v = adj[u][it[u]++];
                if (v == up[0][u] && u != root) continue;
                up[0][v] = u, depth[v] = depth[u] + 1;
                first[v] = tour.size(), tour.push_back(v);
                stack.push_back(v);
                continue;
            }
            last[u] = tour.size(), tour.push_back(u);
            stack.pop_back();
        }

        assert((int)tour.size() == 2 * n);
        for (int k = 1; k <= lg; k++){
            for (int v = 0; v < n; v++) up[k][v] = up[k - 1][up[k - 1][v]];
        }
    }

    int lca(int u, int v) const{
        if (depth[u] < depth[v]) swap(u, v);
        for (int k = lg; k >= 0; k--){
            if (depth[u] - (1 << k) >= depth[v]) u = up[k][u];
        }

        if (u == v) return u;
        for (int k = lg; k >= 0; k--){
            if (up[k][u] != up[k][v]) u = up[k][u], v = up[k][v];
        }
        return up[0][u];
    }

    template <typename Add, typename Remove, typename Get>
    auto run(const vector<pair<int, int>>& paths, Add add, Remove remove, Get get, bool edges = false){
        /// A node inside the tour range twice is not on the path, so each tour position toggles its node
        vector<pair<int, int>> ranges;
        vector<int> extra;
        for (auto [u, v] : paths){
            if (first[u] > first[v]) swap(u, v);
            int l = lca(u, v);
            if (l == u) ranges.push_back({first[u] + edges, first[v]}), extra.push_back(-1);  /// u == v in edge mode is empty
            else ranges.push_back({last[u], first[v]}), extra.push_back(edges ? -1 : l);
        }

        vector<char> inside(n, 0);
        auto toggle = [&](int node){
            if (inside[node]) remove(node);
            else add(node);
            inside[node] ^= 1;
        };

        return mo(2 * n, ranges, [&](int i){ toggle(tour[i]); }, [&](int i){ toggle(tour[i]); }, [&](int i){
            if (extra[i] != -1) toggle(extra[i]);
            auto res = get(i);
            if (extra[i] != -1) toggle(extra[i]);
            return res;
        });
    }
};

int main(){
    /***
     *          0
     *        /   \
     *       1     2
     *      / \     \
     *     3   4     5
    ***/
    TreePathMo tree(6);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}}) tree.add_edge(u, v);
    tree.build(0);
    vector<int> value = {10, 20, 30, 40, 50, 60};
    vector<pair<int, int>> paths = {{3, 4}, {3, 5}, {2, 2}};
    long long sum = 0;
    auto add = [&](int v){ sum += value[v]; };
    auto remove = [&](int v){ sum -= value[v]; };

    auto nodes = tree.run(paths, add, remove, [&](int){ return sum; });
    assert((nodes == vector<long long>{40 + 20 + 50, 40 + 20 + 10 + 30 + 60, 30}));

    auto edges = tree.run(paths, add, remove, [&](int){ return sum; }, true);
    assert((edges == vector<long long>{40 + 50, 40 + 20 + 30 + 60, 0}));  /// edge mode leaves out the lca
    return 0;
}
