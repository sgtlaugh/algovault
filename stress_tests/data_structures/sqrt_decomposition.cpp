#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sqrt_decomposition.cpp"
#undef main

/// Range add and range count below x against a plain array
void check(int n, int ops, long long range){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-range, range);
    SqrtDecomposition<long long> s(a);

    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (stress::rand_int(0, 1)){
            long long v = stress::rand_int(-range, range);
            s.add(l, r, v);
            for (int i = l; i <= r; i++) a[i] += v;
        }
        else{
            long long x = stress::rand_int(0, 3) ? stress::rand_int(-3 * range, 3 * range) : a[stress::rand_int(l, r)];
            int expected = 0;
            for (int i = l; i <= r; i++) expected += a[i] < x;
            assert(s.count_less(l, r, x) == expected);
        }

        int i = stress::rand_int(0, n - 1);
        assert(s.get(i) == a[i]);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++) check(stress::rand_int(1, 80), 200, it % 2 ? 3 : 1000000000);
    for (int n : {1, 2, 3, 4, 15, 16, 17}) check(n, 3000, 5);
    check(100000, 3000, 1000000000);

    /// Values swinging across the whole allowed half range of int, where x - pending would overflow
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 100);
        vector<int> a(n);
        for (auto& x : a) x = stress::rand_int(-1000000000, 1000000000);
        SqrtDecomposition<int> s(a);

        for (int op = 0; op < 200; op++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            int lo = INT_MIN, hi = INT_MAX;
            for (int i = l; i <= r; i++) lo = max(lo, -1000000000 - a[i]), hi = min(hi, 1000000000 - a[i]);
            if (stress::rand_int(0, 1) && lo <= hi){
                int v = stress::rand_int(lo, hi);
                s.add(l, r, v);
                for (int i = l; i <= r; i++) a[i] += v;
            }

            int x = stress::rand_int(-1000000000, 1000000000), expected = 0;
            for (int i = l; i <= r; i++) expected += a[i] < x;
            assert(s.count_less(l, r, x) == expected);
        }
    }

    return 0;
}
