#include "../common.h"

#define main library_main
#include "../../code_library/data_structures/merge_sort_tree.cpp"
#undef main

/// Counts, frequencies and order statistics of both structures against sorting the range directly
void check(const vector<long long>& a, int queries){
    int n = a.size();
    MergeSortTree<long long> tree(a);
    WaveletMatrix<long long> matrix(a);
    long long lo = *min_element(a.begin(), a.end()), hi = *max_element(a.begin(), a.end());

    for (int q = 0; q < queries; q++){
        int l = stress::rand_int(0, n - 1), r = stress::rand_int(l, n - 1);
        if (q % 10 == 0) l = 0, r = n - 1;
        vector<long long> part(a.begin() + l, a.begin() + r + 1);
        sort(part.begin(), part.end());

        long long x = part[stress::rand_int(0, part.size() - 1)];
        if (stress::rand_int(0, 3)) x = stress::rand_int(lo == LLONG_MIN ? lo : lo - 1, hi == LLONG_MAX ? hi : hi + 1);
        int less = lower_bound(part.begin(), part.end(), x) - part.begin();
        int less_equal = upper_bound(part.begin(), part.end(), x) - part.begin();
        assert(tree.count_less(l, r, x) == less && matrix.count_less(l, r, x) == less);
        assert(tree.count_less_equal(l, r, x) == less_equal && matrix.count_less_equal(l, r, x) == less_equal);
        assert(matrix.count_equal(l, r, x) == less_equal - less);

        int k = stress::rand_int(0, r - l);
        assert(tree.kth(l, r, k) == part[k] && matrix.kth(l, r, k) == part[k]);
    }
}

vector<long long> random_array(int n, long long lo, long long hi){
    vector<long long> a(n);
    for (auto& x : a) x = stress::rand_int(lo, hi);
    return a;
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 70);
        long long range = it % 3 == 0 ? 1000000000000000000LL : it % 3 == 1 ? 5 : stress::rand_int(0, 2);
        check(random_array(n, -range, range), 60);
    }

    for (int bits = 0; bits <= 6; bits++){
        vector<long long> a(1 << bits);
        iota(a.begin(), a.end(), 0);
        shuffle(a.begin(), a.end(), stress::rng());
        check(a, 500);
        a.push_back(1 << bits);
        check(a, 500);
    }

    check({LLONG_MIN, LLONG_MAX, LLONG_MIN, 0, LLONG_MAX}, 500);
    for (int n : {1, 2, 3, 127, 128, 129}) check(random_array(n, -3, 3), 2000);
    check(random_array(100000, -1000000000, 1000000000), 300);
    check(random_array(100000, 0, 0), 300);

    return 0;
}
