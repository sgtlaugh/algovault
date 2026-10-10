// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/jump_on_tree
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
    LCA tree(n);
    for (int i = 0; i + 1 < n; i++){
        int a = read_int(), b = read_int();
        tree.add_edge(a, b);
    }
    tree.build(0);

    while (q--){
        int s = read_int(), t = read_int(), i = read_int();
        printf("%d\n", tree.jump(s, t, i));
    }
    return 0;
}
