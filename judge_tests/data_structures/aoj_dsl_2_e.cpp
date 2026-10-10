// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_2_E
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    auto fen = FenwickRangeUpdate<long long>(n);
    while (q--){
        int t;
        if (scanf("%d", &t) != 1) return 0;

        if (t == 0){
            int s, e;
            long long x;
            if (scanf("%d %d %lld", &s, &e, &x) != 3) return 0;
            fen.update(s, e, x);
        }
        else{
            int i;
            if (scanf("%d", &i) != 1) return 0;
            printf("%lld\n", fen.query(i));
        }
    }
    return 0;
}
