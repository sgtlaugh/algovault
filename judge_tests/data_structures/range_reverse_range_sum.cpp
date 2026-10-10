// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_reverse_range_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/treap.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    vector<long long> a(n);
    for (long long& x: a){
        if (scanf("%lld", &x) != 1) return 0;
    }

    ImplicitTreap treap(a);
    while (q--){
        int type, l, r;
        if (scanf("%d %d %d", &type, &l, &r) != 3) return 0;
        if (type == 0 && l < r) treap.reverse(l, r - 1);
        if (type == 1) printf("%lld\n", l < r ? treap.sum(l, r - 1) : 0LL);
    }
    return 0;
}
