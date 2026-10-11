// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/vertex_add_path_sum
// competitive-verifier: TLE 5
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
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        hld.add_edge(u, v);
    }
    hld.build(0);

    FenwickPointUpdate<long long> fen(n);
    for (int v = 0; v < n; v++) fen.update(hld.pos[v], a[v]);

    while (q--){
        int t = read_int(), u = read_int(), v = read_int();
        if (t == 0){
            fen.update(hld.pos[u], v);
            continue;
        }
        long long res = 0;
        for (auto [l, r] : hld.path(u, v)) res += fen.query(l, r);
        printf("%lld\n", res);
    }
    return 0;
}
