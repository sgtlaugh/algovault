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
        vector<long long> a(n), b(n), c(n);

        for (int op = 0; op < 200; op++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(0, n - 1);
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
                assert(point.prefix(l - 1) == accumulate(a.begin(), a.begin() + l, 0LL));
                assert(full.prefix(l - 1) == accumulate(c.begin(), c.begin() + l, 0LL));
            }
        }
    }

    /// lower_bound against a linear prefix scan, zeros are frequent so runs of equal prefix sums occur
    for (long long it = 0; it < stress::scaled(5000); it++){
        int n = it < 70 ? it : stress::rand_int(1, it % 10 ? 40 : 1100);
        long long max_v = it % 3 ? 5 : 1000000000;
        FenwickPointUpdate<long long> fen(n);
        vector<long long> a(n);

        for (int op = 0; op < 100; op++){
            if (n > 0 && stress::rand_int(0, 1)){
                int p = stress::rand_int(0, n - 1);
                long long v = stress::rand_int(0, 1) ? stress::rand_int(0, max_v) : -a[p];
                fen.update(p, v);
                a[p] += v;
                continue;
            }

            long long total = accumulate(a.begin(), a.end(), 0LL);
            long long s = stress::rand_int(-2, total + 2);
            if (n > 0 && stress::rand_int(0, 1)){  /// land on or next to a real prefix sum, which separates < from <=
                int p = stress::rand_int(0, n - 1);
                s = accumulate(a.begin(), a.begin() + p + 1, 0LL) + stress::rand_int(-1, 1);
            }

            int expected = n;
            long long prefix = 0;
            for (int p = 0; p < n; p++){
                prefix += a[p];
                if (prefix >= s){
                    expected = p;
                    break;
                }
            }

            assert(fen.lower_bound(s) == expected);
        }
    }

    /// The O(n) vector constructor against point updates into a zero tree, n = 0 included
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = it < 70 ? it : stress::rand_int(1, it % 10 ? 40 : 1100);
        vector<long long> a(n);
        for (auto& x : a) x = stress::rand_int(-1000000000, 1000000000);
        FenwickPointUpdate<long long> built(a), looped(n);
        for (int i = 0; i < n; i++) looped.update(i, a[i]);
        assert(built.tree == looped.tree);
        for (int p = -1; p < n; p++) assert(built.prefix(p) == accumulate(a.begin(), a.begin() + (p + 1), 0LL));
    }

    return 0;
}
