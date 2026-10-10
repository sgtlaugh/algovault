// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/min_of_mod_of_linear
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/floor_sum.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long n, m, a, b;
        if (scanf("%lld %lld %lld %lld", &n, &m, &a, &b) != 4) return 0;
        printf("%lld\n", min_mod_linear(n, m, a, b));
    }
    return 0;
}
