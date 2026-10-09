#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

/// The default range add, range sum operations against a plain array, through both constructors
int main(){
    for (long long it = 0; it < stress::scaled(15000); it++){
        int n = stress::rand_int(1, it % 20 == 0 ? 2000 : 40);
        vector<long long> a(n, 0);  /// the size constructor starts every element at the identity, 0 for sum
        if (it % 2) for (auto& x : a) x = stress::rand_int(-1000000, 1000000);
        SegmentTree<long long> st = it % 2 ? SegmentTree<long long>(a) : SegmentTree<long long>(n);

        for (int op = 0; op < 100; op++){
            int l = stress::rand_int(1, n), r = stress::rand_int(l, n);
            if (stress::rand_int(0, 1)){
                long long v = stress::rand_int(-1000000, 1000000) * stress::rand_int(0, 1);  /// half the updates add 0
                st.update(l, r, v);
                for (int i = l; i <= r; i++) a[i - 1] += v;
            }
            else assert(st.query(l, r) == accumulate(a.begin() + l - 1, a.begin() + r, 0LL));
        }
    }

    return 0;
}
