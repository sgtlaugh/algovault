#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/mo.cpp"
#undef main

/// Distinct values and sum of every range against direct scans
void check_array(int n, int q, int max_v, int verify = 1 << 30){
    vector<int> a(n);
    for (auto& x : a) x = stress::rand_int(0, max_v);
    vector<pair<int, int>> queries(q);
    for (auto& [l, r] : queries) l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);

    vector<int> freq(max_v + 1, 0);
    int distinct = 0, window = 0;
    long long sum = 0;
    vector<int> got_distinct(q);
    vector<long long> got_sum(q);
    for (int round = 0; round < 2; round++){
        mo(n, queries,
            [&](int i){ distinct += freq[a[i]]++ == 0, sum += a[i], window++; },
            [&](int i){ distinct -= --freq[a[i]] == 0, sum -= a[i], window--; },
            [&](int qi){ got_distinct[qi] = distinct, got_sum[qi] = sum; });
        assert(window == 0 && sum == 0 && distinct == 0);
        for (int i = 0; i < q && i < verify; i++){
            auto [l, r] = queries[i];
            assert(got_sum[i] == accumulate(a.begin() + l, a.begin() + r + 1, 0LL));
            assert(got_distinct[i] == (int)set<int>(a.begin() + l, a.begin() + r + 1).size());
        }
    }
}

vector<pair<int, int>> make_tree(int n, int shape){
    vector<int> label(n);
    iota(label.begin(), label.end(), 0);
    shuffle(label.begin(), label.end(), stress::rng());
    vector<pair<int, int>> edges;
    for (int i = 1; i < n; i++){
        int p = shape == 0 ? stress::rand_int(0, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 0 : (i - 1) / 2;
        edges.push_back({label[i], label[p]});
    }
    return edges;
}

/// Multiset of values on each tree path, nodes or edges, against climbing parents
void check_tree(int n, int q, int shape, int verify = 1 << 30){
    auto edges = make_tree(n, shape);
    int root = stress::rand_int(0, n - 1);
    TreePathMo tree(n);
    for (auto [u, v] : edges) tree.add_edge(u, v);
    tree.build(root);

    vector<int> value(n);
    for (auto& x : value) x = stress::rand_int(0, 20);
    vector<pair<int, int>> paths(q);
    for (auto& [u, v] : paths) u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);

    for (bool edge_mode : {false, true}){
        vector<int> freq(21, 0);
        int distinct = 0;
        long long sum = 0;
        vector<int> got_distinct(q);
        vector<long long> got_sum(q);
        tree.run(paths,
            [&](int v){ distinct += freq[value[v]]++ == 0, sum += value[v]; },
            [&](int v){ distinct -= --freq[value[v]] == 0, sum -= value[v]; },
            [&](int qi){ got_distinct[qi] = distinct, got_sum[qi] = sum; }, edge_mode);
        assert(sum == 0 && distinct == 0);

        for (int i = 0; i < q && i < verify; i++){
            int u = paths[i].first, v = paths[i].second, l;
            for (int x = u, y = v; ; ){
                if (x == y){
                    l = x;
                    break;
                }
                if (tree.depth[x] >= tree.depth[y]) x = tree.up[0][x];
                else y = tree.up[0][y];
            }
            vector<int> nodes;
            for (int x = u; x != l; x = tree.up[0][x]) nodes.push_back(x);
            for (int x = v; x != l; x = tree.up[0][x]) nodes.push_back(x);
            if (!edge_mode) nodes.push_back(l);
            long long s = 0;
            set<int> vals;
            for (int x : nodes) s += value[x], vals.insert(value[x]);
            assert(got_sum[i] == s && got_distinct[i] == (int)vals.size());
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_array(stress::rand_int(1, 60), stress::rand_int(1, 60), it % 2 ? 3 : 50);
        check_tree(stress::rand_int(1, 40), stress::rand_int(1, 40), it % 4);
    }
    check_array(100000, 100000, 1000, 300);
    check_tree(100000, 2000, 1, 200);  /// a path, deep enough to break a recursive Euler tour on an 8 MB stack
    check_tree(30000, 30000, 0, 2000);
    return 0;
}
