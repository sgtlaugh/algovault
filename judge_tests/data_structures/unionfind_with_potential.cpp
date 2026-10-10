// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/unionfind_with_potential
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/disjoint_set.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    auto dsu = WeightedDSU(n, 998244353);
    while (q--){
        int t, u, v;
        if (scanf("%d %d %d", &t, &u, &v) != 3) return 0;
        if (t == 0){
            long long x;
            if (scanf("%lld", &x) != 1) return 0;
            printf("%d\n", dsu.connect(v, u, x) ? 1 : 0);
        }
        else{
            printf("%lld\n", dsu.is_connected(v, u) ? dsu.diff(v, u) : -1LL);
        }
    }
    return 0;
}
