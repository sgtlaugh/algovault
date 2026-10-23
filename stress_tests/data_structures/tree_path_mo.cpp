#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/tree_path_mo.cpp"
#undef main

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

    for (auto mode : {TreePathMo::NODES, TreePathMo::EDGES}){
        vector<int> freq(21, 0);
        int distinct = 0;
        long long sum = 0;
        auto got = tree.run(paths,
            [&](int v){ distinct += freq[value[v]]++ == 0, sum += value[v]; },
            [&](int v){ distinct -= --freq[value[v]] == 0, sum -= value[v]; },
            [&](int qi){ return make_tuple(distinct, sum, qi); }, mode);
        assert(sum == 0 && distinct == 0 && (int)got.size() == q);

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
            if (mode == TreePathMo::NODES) nodes.push_back(l);

            long long s = 0;
            set<int> vals;
            for (int x : nodes) s += value[x], vals.insert(value[x]);
            auto [got_distinct, got_sum, got_qi] = got[i];
            assert(got_qi == i && got_sum == s && got_distinct == (int)vals.size());
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        check_tree(stress::rand_int(1, 40), stress::rand_int(1, 40), it % 4);
    }

    check_tree(100000, 2000, 1, 200);  /// a path, deep enough to break a recursive Euler tour on an 8 MB stack
    check_tree(30000, 30000, 0, 2000);
    return 0;
}
