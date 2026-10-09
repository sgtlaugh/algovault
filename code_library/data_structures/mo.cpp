/***
 *
 * Mo's Algorithm
 * Answers offline range queries by moving a window, for arrays, arrays with point updates, add-only (rollback) state and tree paths
 *
 * Complexity: mo and TreePathMo: O((n + q) * sqrt(n)) calls to add / remove, plus O(q log q) sorting
 *             mo_with_updates: O((q + u) * n^(2/3) + n^(4/3)) calls to add / remove / toggle for u updates, plus O(q log q) sorting
 *             mo_with_rollback: O(n * sqrt(q) + q) calls to add, plus O(q log q) sorting, a rollback costs as much as undoing its adds
 *
 * mo(n, queries, add, remove, answer):
 *   queries[i] = {l, r}, 0-based inclusive ranges over positions 0..n-1
 *   add(i) / remove(i) move position i into / out of the window, answer(qi) records query qi when the window equals it
 *   The window starts and ends empty, so the same callbacks and state can serve several runs
 *
 * mo_with_updates(n, queries, update_pos, add, remove, toggle, answer):
 *   queries[i] = {l, r, t}, the range l..r after the first t updates, 0 <= t <= u
 *   update j sets position update_pos[j], toggle(j) swaps its new value with the array's: alternate calls apply / undo it
 *   Updates are applied in increasing and undone in decreasing order, so the swap is all toggle needs:
 *     toggle = [&](int j){ swap(a[update_pos[j]], value[j]); }
 *   toggle on a position inside the window is wrapped in remove / add, the callbacks only ever read the current array
 *   The run starts and ends with the window empty and no update applied, so the array is left as it was
 *   The n^(2/3) block assumes q and u close to n, with no updates plain mo is faster: n * sqrt(q) against q * n^(2/3)
 *
 * mo_with_rollback(n, queries, add, snapshot, rollback, answer):
 *   For state that supports add but not remove (max, DSU with rollback): queries[i] = {l, r} as in mo
 *   snapshot() saves the state on a stack, rollback() restores the top saved state and pops it, at most 2 are saved
 *   The window is always contiguous and add(i) extends it at one end, it starts and ends empty
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

template <typename Add, typename Remove, typename Answer>
void mo(int n, const vector<pair<int, int>>& queries, Add add, Remove remove, Answer answer){
    int q = queries.size();
    if (q == 0) return;
    int block = max(1, (int)(n / sqrt((double)q)));

    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? queries[a].second > queries[b].second : queries[a].second < queries[b].second;
    });

    int cur_l = 0, cur_r = -1;
    for (int qi : order){
        auto [l, r] = queries[qi];
        assert(0 <= l && l <= r && r < n);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        answer(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);  /// leave the window empty so the caller's state can be reused
}

template <typename Add, typename Remove, typename Toggle, typename Answer>
void mo_with_updates(int n, const vector<array<int, 3>>& queries, const vector<int>& update_pos, Add add, Remove remove, Toggle toggle, Answer answer){
    int q = queries.size(), u = update_pos.size();
    if (q == 0) return;
    int block = max(1, (int)cbrt((double)n * n));

    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int la = queries[a][0] / block, lb = queries[b][0] / block;
        if (la != lb) return la < lb;
        int ra = queries[a][1] / block, rb = queries[b][1] / block;
        if (ra != rb) return ra < rb;
        return queries[a][2] < queries[b][2];
    });

    int cur_l = 0, cur_r = -1, cur_t = 0;
    auto flip = [&](int j){
        int p = update_pos[j];
        assert(0 <= p && p < n);
        bool inside = cur_l <= p && p <= cur_r;
        if (inside) remove(p);
        toggle(j);
        if (inside) add(p);
    };

    for (int qi : order){
        auto [l, r, t] = queries[qi];
        assert(0 <= l && l <= r && r < n && 0 <= t && t <= u);
        while (cur_t < t) flip(cur_t++);
        while (cur_t > t) flip(--cur_t);
        while (cur_l > l) add(--cur_l);
        while (cur_r < r) add(++cur_r);
        while (cur_l < l) remove(cur_l++);
        while (cur_r > r) remove(cur_r--);
        answer(qi);
    }

    while (cur_r >= cur_l) remove(cur_r--);
    while (cur_t > 0) toggle(--cur_t);
}

template <typename Add, typename Snapshot, typename Rollback, typename Answer>
void mo_with_rollback(int n, const vector<pair<int, int>>& queries, Add add, Snapshot snapshot, Rollback rollback, Answer answer){
    int q = queries.size();
    if (q == 0) return;
    int block = max(1, (int)(n / sqrt((double)q)));

    vector<int> order(q);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        int ba = queries[a].first / block, bb = queries[b].first / block;
        if (ba != bb) return ba < bb;
        return queries[a].second < queries[b].second;
    });

    /// Per block of l the right part grows from the block end and is kept, the left part is added and rolled back per query
    for (int i = 0; i < q; ){
        int b = queries[order[i]].first / block, block_end = min(n - 1, (b + 1) * block - 1), cur_r = block_end;
        snapshot();
        for (; i < q && queries[order[i]].first / block == b; i++){
            int qi = order[i];
            auto [l, r] = queries[qi];
            assert(0 <= l && l <= r && r < n);
            while (cur_r < r) add(++cur_r);
            snapshot();
            if (r <= block_end){
                for (int j = l; j <= r; j++) add(j);  /// sorted by r, so these run before the right part grows
            }
            else{
                for (int j = block_end; j >= l; j--) add(j);
            }
            answer(qi);
            rollback();
        }
        rollback();
    }
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
        });
    }
};

int main(){
    vector<int> a = {1, 2, 1, 3, 2, 2, 4};
    vector<pair<int, int>> queries = {{0, 6}, {0, 2}, {3, 5}, {4, 4}, {2, 3}};
    vector<int> freq(5, 0), distinct(queries.size());
    int count = 0;
    mo((int)a.size(), queries,
        [&](int i){ count += freq[a[i]]++ == 0; },
        [&](int i){ count -= --freq[a[i]] == 0; },
        [&](int qi){ distinct[qi] = count; });
    assert((distinct == vector<int>{4, 2, 2, 1, 2}));

    vector<int> b = {1, 2, 1, 3, 2}, update_pos = {1, 3}, update_value = {1, 2};
    vector<array<int, 3>> timed = {{0, 4, 0}, {0, 4, 1}, {0, 4, 2}, {1, 3, 0}, {1, 3, 2}, {0, 2, 1}, {3, 3, 1}, {3, 3, 2}};
    vector<int> timed_distinct(timed.size());
    mo_with_updates((int)b.size(), timed, update_pos,
        [&](int i){ count += freq[b[i]]++ == 0; },
        [&](int i){ count -= --freq[b[i]] == 0; },
        [&](int j){ swap(b[update_pos[j]], update_value[j]); },
        [&](int qi){ timed_distinct[qi] = count; });
    assert((timed_distinct == vector<int>{3, 3, 2, 3, 2, 1, 1, 1}));
    assert((b == vector<int>{1, 2, 1, 3, 2} && update_value == vector<int>{1, 2} && count == 0));

    vector<int> c = {9, 8, 7, 8, 9, 4, 8}, c_freq(10, 0);
    vector<pair<int, int>> ranges = {{0, 6}, {1, 6}, {1, 3}, {0, 4}, {5, 5}, {2, 4}};
    vector<pair<int, long long>> history;
    vector<size_t> marks;
    vector<long long> importance(ranges.size());
    long long best = 0;
    mo_with_rollback((int)c.size(), ranges,
        [&](int i){
            history.push_back({c[i], best});
            best = max(best, (long long)c[i] * ++c_freq[c[i]]);
        },
        [&](){ marks.push_back(history.size()); },
        [&](){
            for (; history.size() > marks.back(); history.pop_back()) c_freq[history.back().first]--, best = history.back().second;
            marks.pop_back();
        },
        [&](int qi){ importance[qi] = best; });
    assert((importance == vector<long long>{24, 24, 16, 18, 4, 9}));
    assert(best == 0 && history.empty() && marks.empty());

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
    vector<pair<int, int>> paths = {{3, 4}, {3, 5}, {1, 3}, {2, 2}, {5, 0}};

    vector<long long> sums(paths.size());
    long long sum = 0;
    tree.run(paths, [&](int v){ sum += value[v]; }, [&](int v){ sum -= value[v]; }, [&](int qi){ sums[qi] = sum; });
    assert((sums == vector<long long>{40 + 20 + 50, 40 + 20 + 10 + 30 + 60, 20 + 40, 30, 60 + 30 + 10}));

    tree.run(paths, [&](int v){ sum += value[v]; }, [&](int v){ sum -= value[v]; }, [&](int qi){ sums[qi] = sum; }, true);
    assert((sums == vector<long long>{40 + 50, 40 + 20 + 30 + 60, 40, 0, 60 + 30}));

    return 0;
}
