#include "common.h"

#define main library_main
#include "../code_library/merge_sort_tree.cpp"
#undef main

/// Counts and order statistics against sorting the range directly
void check(int n, int queries, long long range){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(-range, range);
    MergeSortTree<long long> tree(a);
    for (int q = 0; q < queries; q++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (q % 10 == 0) l = 0, r = n - 1;
        vector<long long> part(a.begin() + l, a.begin() + r + 1);
        sort(part.begin(), part.end());
        long long x = stress::rand_int(0, 3) ? stress::rand_int(-range - 1, range + 1) : part[stress::rand_int(0, part.size() - 1)];
        assert(tree.count_less(l, r, x) == lower_bound(part.begin(), part.end(), x) - part.begin());
        assert(tree.count_less_equal(l, r, x) == upper_bound(part.begin(), part.end(), x) - part.begin());
        int k = stress::rand_int(0, r - l);
        assert(tree.kth(l, r, k) == part[k]);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        check(stress::rand_int(1, 70), 60, it % 3 ? 5 : 1000000000000000000LL);
    }
    for (int n : {1, 2, 3, 127, 128, 129}) check(n, 2000, 3);
    check(100000, 300, 1000000000);
    return 0;
}
