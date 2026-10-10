// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_5_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/dynamic_diameter.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    DynamicDiameter tree(n);
    for (int i = 0; i + 1 < n; i++){
        int s, t, w;
        if (scanf("%d %d %d", &s, &t, &w) != 3) return 0;
        tree.add_edge(s, t, w);
    }
    tree.build();

    printf("%lld\n", tree.diameter());
    return 0;
}
