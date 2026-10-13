#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/disjoint_set.cpp"
#undef main

/// Union by size keeps every node within log2(component size) parent hops of its root
template <typename Dsu>
void check_depths(const Dsu& dsu){
    for (int i = 0; i < (int)dsu.parent.size(); i++){
        int depth = 0, root = i;
        while (dsu.parent[root] != root) root = dsu.parent[root], depth++;
        assert(depth < 31 && (1 << depth) <= dsu.sz[root]);
    }
}

/// Applies one union to a label array in O(n)
void relabel(vector<int>& label, int a, int b){
    int from = label[a], to = label[b];
    if (from != to) for (auto& x : label) if (x == from) x = to;
}

/// Component labels recomputed from scratch over the edge list
vector<int> naive_labels(int n, const vector<pair<int, int>>& edges){
    vector<int> label(n + 1);
    iota(label.begin(), label.end(), 0);
    for (auto [a, b] : edges) relabel(label, a, b);
    return label;
}

/// Random unions and queries against a naive label array, merged by relabeling
int main(){
    for (int n : {1, 2, 1000, 100000}){
        DSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));  /// links the grown root under a fresh node each time
        check_depths(chain);
        assert(chain.size(0) == n + 1 && !chain.connect(0, n));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 400);
        DSU dsu(n);
        vector<int> label = naive_labels(n, {});

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            if (stress::rand_int(0, 2) == 0){
                assert(dsu.connect(a, b) == (label[a] != label[b]));
                relabel(label, a, b);
            }
            assert(dsu.same(a, b) == (label[a] == label[b]));
            assert(dsu.size(a) == count(label.begin(), label.end(), label[a]));
        }
        check_depths(dsu);
    }

    return 0;
}
