#include "common.h"

#define main library_main
#include "../code_library/lca.cpp"
#undef main

struct NaiveTree{
    vector<int> parent, level;
    vector<long long> up_weight;

    NaiveTree(int n, int root, const vector<array<long long, 3>>& edges) : parent(n, -1), level(n, 0), up_weight(n, 0){
        vector<vector<pair<int, long long>>> adj(n);
        for (auto [u, v, w] : edges) adj[u].push_back({(int)v, w}), adj[v].push_back({(int)u, w});
        vector<int> queue = {root};
        vector<bool> seen(n, false);
        seen[root] = true;
        for (int i = 0; i < (int)queue.size(); i++){
            int u = queue[i];
            for (auto [v, w] : adj[u]){
                if (seen[v]) continue;
                seen[v] = true, parent[v] = u, level[v] = level[u] + 1, up_weight[v] = w;
                queue.push_back(v);
            }
        }
    }

    int lca(int u, int v) const{
        while (level[u] > level[v]) u = parent[u];
        while (level[v] > level[u]) v = parent[v];
        while (u != v) u = parent[u], v = parent[v];
        return u;
    }

    long long dist(int u, int v) const{
        long long res = 0;
        while (level[u] > level[v]) res += up_weight[u], u = parent[u];
        while (level[v] > level[u]) res += up_weight[v], v = parent[v];
        while (u != v) res += up_weight[u] + up_weight[v], u = parent[u], v = parent[v];
        return res;
    }

    int kth(int v, int k) const{
        if (k < 0 || k > level[v]) return -1;
        while (k--) v = parent[v];
        return v;
    }
};

/// Random tree of a given shape with random labels and weights
vector<array<long long, 3>> make_tree(int n, int shape, long long max_w){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());

    vector<array<long long, 3>> edges;
    for (int i = 1; i < n; i++){
        int p;
        if (shape == 0) p = stress::rand_int(0, i - 1);
        else if (shape == 1) p = i - 1;
        else if (shape == 2) p = 0;
        else if (shape == 3) p = i % 2 ? i - 1 : max(0, i - 2);
        else p = i < n / 2 ? i - 1 : stress::rand_int(max(0, n / 2 - 3), i - 1);
        edges.push_back({label[i], label[p], stress::rand_int(0, max_w)});
    }
    shuffle(edges.begin(), edges.end(), stress::rng());
    return edges;
}

void check(int n, int shape, int queries, long long max_w){
    auto edges = make_tree(n, shape, max_w);
    int root = stress::rand_int(0, n - 1);

    LCA a(n);
    LinearLCA b(n);
    for (auto [u, v, w] : edges) a.add_edge(u, v, w), b.add_edge(u, v, w);
    a.build(root), b.build(root);
    NaiveTree naive(n, root, edges);

    for (int q = 0; q < queries; q++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        int l = naive.lca(u, v);
        assert(a.lca(u, v) == l);
        assert(b.lca(u, v) == l);
        long long d = naive.dist(u, v);
        assert(a.dist(u, v) == d && b.dist(u, v) == d);
        assert(a.depth(u) == naive.level[u] && b.depth(u) == naive.level[u]);
        int k = stress::rand_int(-1, naive.level[u] + 1);
        assert(a.kth_ancestor(u, k) == naive.kth(u, k));
    }
}

/// LinearRMQ against a direct scan, sizes straddling the 64 element blocks and many equal values
void check_rmq(int n, int max_v){
    vector<int> val(n);
    for (auto& x : val) x = stress::rand_int(-max_v, max_v);
    LinearRMQ rmq(val);
    for (int q = 0; q < 300; q++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (q < 64) r = min(n - 1, l + q);
        assert(val[rmq.query(l, r)] == *min_element(val.begin() + l, val.begin() + r + 1));
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 200 ? it % 40 + 1 : stress::rand_int(1, 300);
        check(n, it % 5, 200, it % 3 ? 1000000000LL : 1);
    }

    for (long long it = 0; it < stress::scaled(10); it++){
        check(stress::rand_int(20000, 60000), it % 5, 3000, 1000000000LL);
    }
    check(200000, 1, 100, 1000000000LL);  /// a path, deep enough to overflow an 8 MB stack with recursive traversals

    for (int n : {1, 2, 63, 64, 65, 127, 128, 129, 1000}){
        for (int it = 0; it < 20; it++) check_rmq(n, it % 2 ? 3 : 1000000000);
    }
    for (long long it = 0; it < stress::scaled(200); it++) check_rmq(stress::rand_int(1, 5000), stress::rand_int(0, 1) ? 5 : 1000000000);

    return 0;
}
