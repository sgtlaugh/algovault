/***
 *
 * Mo's Algorithm On Tree Paths
 * Offline path queries on a tree, by Mo's algorithm over an Euler tour that lists each node twice
 *
 * Complexity:
 *   - build: O(n log n) time and memory
 *   - lca: O(log n)
 *   - run: O(n sqrt(q) + q) calls to add and remove, plus O(q log q) sorting and O(q log n) for the lcas
 *
 * TreePathMo tree(n): nodes 0..n - 1 and no edges yet
 * tree.add_edge(u, v): adds the tree edge between u and v
 * tree.build(root = 0): roots the tree at root once all n - 1 edges are added
 * tree.lca(u, v): the lowest common ancestor of u and v, after build
 * tree.run(paths, add, remove, get, mode = NODES): the answers, res[i] = get(i) while the window holds the path paths[i] = {u, v}
 *
 * add(node) and remove(node) receive tree nodes, get(i) returns the answer for the current window
 * mode = NODES: the window holds every node on the path u..v
 * mode = EDGES: the window holds the path without the lca, so store the value of edge (parent[v], v) on v
 * The window starts and ends empty, so the same callbacks and state can serve several runs
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

    while (cur_r >= cur_l) remove(cur_r--);
    return res;
}
/// END COPY mo

struct TreePathMo{
    enum Mode{NODES, EDGES};

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
                if (v == up[0][u]) continue;
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
    auto run(const vector<pair<int, int>>& paths, Add add, Remove remove, Get get, Mode mode = NODES){
        /// A node inside the tour range twice is not on the path, so each tour position toggles its node
        /// With first[u] <= first[v] and u an ancestor of v, [first[u], first[v]] holds exactly the path u..v
        /// Otherwise [last[u], first[v]] holds the path minus the lca, whose entry and exit both lie outside it
        vector<pair<int, int>> ranges;
        vector<int> lca_to_add;  /// toggled in around get(i) for node paths whose range misses the lca, else -1
        for (auto [u, v] : paths){
            if (first[u] > first[v]) swap(u, v);
            int w = lca(u, v);
            if (w == u){
                ranges.push_back({first[u] + (mode == EDGES), first[v]});  /// edge mode drops the lca u, u == v leaves it empty
                lca_to_add.push_back(-1);
            }
            else{
                ranges.push_back({last[u], first[v]});
                lca_to_add.push_back(mode == EDGES ? -1 : w);
            }
        }

        vector<char> inside(n, 0);
        auto toggle = [&](int node){
            if (inside[node]) remove(node);
            else add(node);
            inside[node] ^= 1;
        };

        return mo(2 * n, ranges, [&](int i){ toggle(tour[i]); }, [&](int i){ toggle(tour[i]); }, [&](int i){
            if (lca_to_add[i] != -1) toggle(lca_to_add[i]);
            auto res = get(i);
            if (lca_to_add[i] != -1) toggle(lca_to_add[i]);
            return res;
        });
    }
};

int main(){
    /// Path sums over nodes, then over edges with each edge's value stored on its lower node
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

    auto edges = tree.run(paths, add, remove, [&](int){ return sum; }, TreePathMo::EDGES);
    assert((edges == vector<long long>{40 + 50, 40 + 20 + 30 + 60, 0}));  /// edge mode leaves out the lca

    return 0;
}
