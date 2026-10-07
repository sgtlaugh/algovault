#include "common.h"

#define main library_main
#include "../code_library/disjoint_set.cpp"
#undef main

/// Random unions and queries against a naive label array, merged by relabeling
int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(0, it % 10 ? 15 : 400);
        DSU dsu(n);
        vector<int> label(n + 1);
        iota(label.begin(), label.end(), 0);

        for (int q = stress::rand_int(0, 4 * n + 5); q; q--){
            int a = stress::rand_int(0, n), b = stress::rand_int(0, n);
            if (stress::rand_int(0, 2) == 0){
                dsu.connect(a, b);
                int from = label[a], to = label[b];
                for (auto& x : label) if (x == from) x = to;
            }
            assert(dsu.is_connected(a, b) == (label[a] == label[b]));
            assert(dsu.component_size(a) == count(label.begin(), label.end(), label[a]));
        }
    }
    return 0;
}
