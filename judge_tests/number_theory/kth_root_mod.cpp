// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/kth_root_mod
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/modular_roots.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long k, y, p;
        if (scanf("%lld %lld %lld", &k, &y, &p) != 3) return 0;
        printf("%lld\n", kth_root(y, k, p));
    }
    return 0;
}
