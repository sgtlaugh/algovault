// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/primitive_root
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/modular_roots.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;
    while (q--){
        long long p;
        if (scanf("%lld", &p) != 1) return 0;
        printf("%lld\n", primitive_root(p));
    }
    return 0;
}
