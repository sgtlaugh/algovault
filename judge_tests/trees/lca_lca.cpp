// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lca
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/lca.cpp"
#undef main

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    LinearLCA tree(n);
    for (int i = 1; i < n; i++) tree.add_edge(read_int(), i);
    tree.build(0);

    while (q--){
        int u = read_int(), v = read_int();
        printf("%d\n", tree.lca(u, v));
    }
    return 0;
}
