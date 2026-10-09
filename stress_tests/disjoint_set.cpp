#include "common.h"

#define main library_main
#include "../code_library/disjoint_set.cpp"
#undef main

/// Union by size keeps every node within log2(component size) parent hops of its root
void check_depths(const DSU& dsu){
    for (int i = 0; i < (int)dsu.parent.size(); i++){
        int depth = 0, root = i;
        while (dsu.parent[root] != root) root = dsu.parent[root], depth++;
        assert((1LL << depth) <= dsu.counter[root]);
    }
}

/// Random unions and queries against a naive label array, merged by relabeling
int main(){
    for (int n : {1, 2, 1000, 100000}){
        DSU chain(n);
        for (int i = 0; i < n; i++) assert(chain.connect(i, i + 1));  /// links the grown root under a fresh node each time
        check_depths(chain);
        assert(chain.component_size(0) == n + 1 && !chain.connect(0, n));
    }

    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 400);
        DSU dsu(n);
        vector<int> label(n + 1);
        iota(label.begin(), label.end(), 0);

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            if (stress::rand_int(0, 2) == 0){
                int from = label[a], to = label[b];
                assert(dsu.connect(a, b) == (from != to));
                for (auto& x : label) if (x == from) x = to;
            }
            assert(dsu.is_connected(a, b) == (label[a] == label[b]));
            assert(dsu.component_size(a) == count(label.begin(), label.end(), label[a]));
        }
        check_depths(dsu);
    }
    return 0;
}
