// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_6_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/hlpp.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    HLPP g(n, 0, n - 1);
    for (int i = 0; i < m; i++){
        int u, v;
        long long c;
        if (scanf("%d %d %lld", &u, &v, &c) != 3) return 0;
        g.add_directed_edge(u, v, c);
    }

    printf("%lld\n", g.maxflow());
    return 0;
}
