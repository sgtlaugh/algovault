// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_5_E
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/hld.cpp"
#undef main
#define main fenwick_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    HLD hld(n);
    for (int u = 0; u < n; u++){
        int k, c;
        if (scanf("%d", &k) != 1) return 0;
        while (k--){
            if (scanf("%d", &c) != 1) return 0;
            hld.add_edge(u, c);
        }
    }
    hld.build(0);

    FenwickFull<long long> fen(n);
    int q;
    if (scanf("%d", &q) != 1) return 0;
    while (q--){
        int t, v, w;
        if (scanf("%d %d", &t, &v) != 2) return 0;
        if (t == 0){
            if (scanf("%d", &w) != 1) return 0;
            hld.update_path(fen, 0, v, (long long)w, true);
            continue;
        }
        long long res = 0;
        for (auto [l, r] : hld.path(0, v, true)) res += fen.query(l, r);
        printf("%lld\n", res);
    }
    return 0;
}
