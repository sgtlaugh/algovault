// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/minimum_steiner_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/steiner_tree.cpp"
#undef main

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    SteinerTree g(n);
    vector<long long> weight(m);
    for (int i = 0; i < m; i++){
        int u, v;
        if (scanf("%d %d %lld", &u, &v, &weight[i]) != 3) return 0;
        g.add_edge(u, v, weight[i]);
    }

    int k;
    if (scanf("%d", &k) != 1) return 0;
    vector<int> terminals(k);
    for (int& t : terminals){
        if (scanf("%d", &t) != 1) return 0;
    }

    vector<int> edges = g.tree(terminals);
    long long total = 0;
    for (int id : edges) total += weight[id];
    assert(total == g.solve(terminals));

    printf("%lld %d\n", total, (int)edges.size());
    for (int id : edges) printf("%d ", id);
    printf("\n");
    return 0;
}
