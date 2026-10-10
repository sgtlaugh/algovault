// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sqrt_mod
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/modular_roots.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long y, p;
        if (scanf("%lld %lld", &y, &p) != 2) return 0;
        printf("%lld\n", sqrt_mod(y, p));
    }
    return 0;
}
