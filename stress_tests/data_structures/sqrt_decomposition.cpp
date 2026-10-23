#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/sqrt_decomposition.cpp"
#undef main

/// Range add and range count below x against a plain array, adds at most max_add_length long when it is positive
void check(int n, int ops, long long range, int max_add_length = 0){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-range, range);
    SqrtDecomposition<long long> s(a);

    for (int op = 0; op < ops; op++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (stress::rand_int(0, 1)){
            if (max_add_length > 0 && stress::rand_int(0, 3)) r = min(n - 1, l + (int)stress::rand_int(0, max_add_length - 1));
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

/// Blocks pushed from one edge of the allowed half range of T to the other, where x - pending would overflow
template<typename T>
void check_extremes(T bound){
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 100);
        vector<T> a(n);
        for (auto& x : a) x = it % 2 ? stress::rand_int(-bound, bound) : -bound + stress::rand_int(0, 3);
        SqrtDecomposition<T> s(a);

        for (int op = 0; op < 200; op++){
            int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
            T lo = numeric_limits<T>::min(), hi = numeric_limits<T>::max();
            for (int i = l; i <= r; i++) lo = max(lo, -bound - a[i]), hi = min(hi, bound - a[i]);
            if (stress::rand_int(0, 1) && lo <= hi){
                int edge = stress::rand_int(0, 2);
                T v = edge == 0 ? lo : edge == 1 ? hi : stress::rand_int(lo, hi);
                s.add(l, r, v);
                for (int i = l; i <= r; i++) a[i] += v;
            }

            T x = stress::rand_int(0, 3) ? stress::rand_int(-bound, bound) : stress::rand_int(0, 1) ? -bound : bound;
            int expected = 0;
            for (int i = l; i <= r; i++) expected += a[i] < x;
            assert(s.count_less(l, r, x) == expected);
        }
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++) check(stress::rand_int(1, 80), 200, it % 2 ? 3 : 1000000000);
    for (int n : {1, 2, 3, 4, 15, 16, 17}) check(n, 3000, 5);
    check(100000, 3000, 1000000000);

    /// Point and short adds, mostly inside one block, interleaved with whole-block adds and full-range counts
    for (long long it = 0; it < stress::scaled(2000); it++) check(stress::rand_int(1, 120), 300, it % 2 ? 4 : 1000000000, stress::rand_int(1, 3));
    check(200000, 5000, 1000000000, 1);

    check_extremes<int>(1000000000);
    check_extremes<long long>(LLONG_MAX / 2);

    return 0;
}
