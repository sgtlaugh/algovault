// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/discrete_logarithm_mod
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/discrete_logarithm.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        int x, y, m;
        if (scanf("%d %d %d", &x, &y, &m) != 3) return 0;
        printf("%d\n", discrete_log(x, y, m));
    }
    return 0;
}
