// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_kth_smallest
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/merge_sort_tree.cpp"
#undef main

/// The same problem verifies RangeKth, so all three structures answer every query
#define main persistent_segment_tree_main
#include "../../code_library/data_structures/persistent_segment_tree.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    vector<int> a(n);
    for (int& x: a){
        if (scanf("%d", &x) != 1) return 0;
    }

    MergeSortTree<int> merge_sort_tree(a);
    WaveletMatrix<int> wavelet(a);
    RangeKth<int> range_kth(a);
    while (q--){
        int l, r, k;
        if (scanf("%d %d %d", &l, &r, &k) != 3) return 0;

        int x = merge_sort_tree.kth(l, r - 1, k), y = wavelet.kth(l, r - 1, k), z = range_kth.kth(l, r - 1, k);
        if (x != y || x != z){
            fprintf(stderr, "kth(%d, %d, %d): MergeSortTree %d, WaveletMatrix %d, RangeKth %d\n", l, r - 1, k, x, y, z);
            return 1;
        }
        printf("%d\n", x);
    }
    return 0;
}
