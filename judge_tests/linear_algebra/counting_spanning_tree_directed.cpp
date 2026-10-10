// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/counting_spanning_tree_directed
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/determinant.cpp"
#undef main

int main(){
    int n, m, r;
    if (scanf("%d %d %d", &n, &m, &r) != 3) return 0;
    Arborescences d(n);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d", &u, &v) != 2) return 0;
        d.add_edge(u, v);
    }

    printf("%lld\n", d.count(r, 998244353));
    return 0;
}
