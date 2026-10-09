#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

/// All three variants against a plain array, interleaving updates and queries
int main(){
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = stress::rand_int(1, it % 10 ? 20 : 300);
        FenwickPointUpdate<long long> point(n);
        FenwickRangeUpdate<long long> range(n);
        FenwickFull<long long> full(n);
        vector<long long> a(n + 1), b(n + 1), c(n + 1);

        for (int op = 0; op < 200; op++){
            int l = stress::rand_int(1, n), r = stress::rand_int(1, n);
            long long v = stress::rand_int(-1000000000, 1000000000);
            if (l > r) swap(l, r);
            if (op % 7 == 0) swap(l, r);  /// l > r is an empty range, updates must be no-ops and queries 0

            if (stress::rand_int(0, 1)){
                point.update(l, v);
                a[l] += v;
                range.update(l, r, v);
                full.update(l, r, v);
                for (int i = l; i <= r; i++) b[i] += v, c[i] += v;
            }
            else{
                assert(point.query(l, r) == accumulate(a.begin() + l, a.begin() + max(l, r + 1), 0LL));
                assert(range.query(l) == b[l]);
                assert(full.query(l, r) == accumulate(c.begin() + l, c.begin() + max(l, r + 1), 0LL));
                assert(full.query(l) == accumulate(c.begin() + 1, c.begin() + l + 1, 0LL));
            }
        }
    }

    return 0;
}
