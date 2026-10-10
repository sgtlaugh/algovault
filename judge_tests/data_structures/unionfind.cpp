// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/unionfind
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/disjoint_set.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    auto dsu = DSU(n);

    while (q--){
        int t = read_int(), u = read_int(), v = read_int();
        if (t == 0) dsu.connect(u, v);
        else putchar(dsu.is_connected(u, v) ? '1' : '0'), putchar('\n');
    }
    return 0;
}
