/***
 *
 * Mo's Algorithm
 * Answers offline range queries by moving a window, for arrays and tree paths
 * mo_with_updates.cpp handles arrays that also receive point assignments, mo_with_rollback.cpp state that cannot remove
 *
 * Complexity: O((n + q) * sqrt(n)) calls to add / remove, plus O(q log q) sorting
 *
 * mo(n, queries, add, remove, get) returns the answers, res[qi] = get(qi) taken while the window is exactly query qi:
 *   queries[i] = {l, r}, 0-based inclusive ranges over positions 0..n-1, r = l - 1 is an empty range
 *   add(i) / remove(i) move position i into / out of the window, get(qi) returns the answer for the current window
 *   The window starts and ends empty, so the same callbacks and state can serve several runs
 *
 *   auto distinct = mo(n, queries, [&](int i){ count += freq[a[i]]++ == 0; },
 *                      [&](int i){ count -= --freq[a[i]] == 0; }, [&](int){ return count; });
 *
 * TreePathMo tree(n); tree.add_edge(u, v); tree.build(root);
 * tree.run(paths, add, remove, answer, edges):
 *   paths[i] = {u, v}, add(node) / remove(node) receive tree nodes
 *   edges = false: the window holds every node on the path u..v
 *   edges = true: the window holds every node on the path except the lca, store the value of edge (parent[v], v) on v
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename Add, typename Remove, typename Get>
auto mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Get get){
    int q = queries.size();
    vector<decay_t<decltype(get(0))>> res(q);
    if (q == 0) return res;
    int block = max(1, (int)(n / sqrt((double)q)));

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

    template <typename Add, typename Remove, typename Answer>
    void run(const vector<pair<int, int>>& paths, Add add, Remove remove, Answer answer, bool edges = false){
        /// A node inside the tour range twice is not on the path, so each tour position toggles its node
        vector<pair<int, int>> ranges;
        vector<int> extra;
        for (auto [u, v] : paths){
            if (first[u] > first[v]) swap(u, v);
            int l = lca(u, v);
            if (l == u) ranges.push_back({first[u] + edges, first[v]}), extra.push_back(-1);
            else ranges.push_back({last[u], first[v]}), extra.push_back(edges ? -1 : l);
        }

        vector<char> inside(n, 0);
        auto toggle = [&](int node){
            if (inside[node]) remove(node);
            else add(node);
            inside[node] ^= 1;
        };

        vector<pair<int, int>> real;
        vector<int> index;
        for (int i = 0; i < (int)ranges.size(); i++){
            if (ranges[i].first <= ranges[i].second) real.push_back(ranges[i]), index.push_back(i);
            else answer(i);  /// u == v in edge mode, the path holds no edges
        }

        mo(2 * n, real, [&](int i){ toggle(tour[i]); }, [&](int i){ toggle(tour[i]); }, [&](int qi){
            int i = index[qi];
            if (extra[i] != -1) toggle(extra[i]);
            answer(i);
            if (extra[i] != -1) toggle(extra[i]);
            return 0;
        });
    }
};

int main(){
    /// Distinct values per range: the callbacks keep a frequency table and a running count
    vector<int> a = {1, 2, 1, 3, 2, 2, 4};
    vector<pair<int, int>> queries = {{0, 6}, {0, 2}, {3, 5}, {4, 3}};
    vector<int> freq(5, 0);
    int count = 0;
    auto distinct = mo((int)a.size(), queries,
        [&](int i){ count += freq[a[i]]++ == 0; },
        [&](int i){ count -= --freq[a[i]] == 0; },
        [&](int){ return count; });
    assert((distinct == vector<int>{4, 2, 2, 0}));  /// {4, 3} is empty

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
    vector<pair<int, int>> paths = {{3, 4}, {3, 5}};
    vector<long long> sums(paths.size());
    long long sum = 0;
    tree.run(paths, [&](int v){ sum += value[v]; }, [&](int v){ sum -= value[v]; }, [&](int qi){ sums[qi] = sum; });
    assert((sums == vector<long long>{40 + 20 + 50, 40 + 20 + 10 + 30 + 60}));

    tree.run(paths, [&](int v){ sum += value[v]; }, [&](int v){ sum -= value[v]; }, [&](int qi){ sums[qi] = sum; }, true);
    assert((sums == vector<long long>{40 + 50, 40 + 20 + 30 + 60}));  /// edge mode leaves out the lca
    return 0;
}
