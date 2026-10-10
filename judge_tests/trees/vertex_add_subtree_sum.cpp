// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/vertex_add_subtree_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/hld.cpp"
#undef main
#define main fenwick_main
#include "../../code_library/data_structures/fenwick_tree.cpp"
#undef main

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<int> a(n);
    for (int& x : a) x = read_int();
    HLD hld(n);
    for (int i = 1; i < n; i++) hld.add_edge(read_int(), i);
    hld.build(0);

    FenwickPointUpdate<long long> fen(n);
    for (int v = 0; v < n; v++) fen.update(hld.pos[v], a[v]);

    while (q--){
        int t = read_int(), u = read_int();
        if (t == 0){
            fen.update(hld.pos[u], read_int());
            continue;
        }
        auto [l, r] = hld.subtree(u);
        printf("%lld\n", fen.query(l, r));
    }
    return 0;
}
