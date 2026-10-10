// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/GRL_1_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/bellman_ford.cpp"
#undef main

int main(){
    int n, m, r;
    if (scanf("%d %d %d", &n, &m, &r) != 3) return 0;
    BellmanFord g(n);
    for (int i = 0; i < m; i++){
        int s, t;
        long long d;
        if (scanf("%d %d %lld", &s, &t, &d) != 3) return 0;
        g.add_edge(s, t, d);
    }

    auto dist = g.shortest_paths(r);
    if (count(dist.begin(), dist.end(), BellmanFord::NEG_INF)){
        puts("NEGATIVE CYCLE");
        return 0;
    }

    for (long long d : dist){
        if (d == BellmanFord::INF) puts("INF");
        else printf("%lld\n", d);
    }
    return 0;
}
