// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_2_G
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main fenwick_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main
#define main segment_tree_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    auto fen = FenwickFull<long long>(n);
    auto seg = SegmentTree<long long>(n);
    while (q--){
        int t, s, e;
        if (scanf("%d %d %d", &t, &s, &e) != 3) return 0;
        if (t == 0){
            long long x;
            if (scanf("%lld", &x) != 1) return 0;
            fen.update(s - 1, e - 1, x);
            seg.update(s - 1, e - 1, x);
        }
        else{
            long long res = fen.query(s - 1, e - 1);
            assert(seg.query(s - 1, e - 1) == res);
            printf("%lld\n", res);
        }
    }
    return 0;
}
