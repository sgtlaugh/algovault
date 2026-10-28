// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/minimum_steiner_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/steiner_tree.cpp"
#undef main

int main(){
    int n, m;
    scanf("%d %d", &n, &m);
    SteinerTree g(n);
    vector<long long> weight(m);
    for (int i = 0; i < m; i++){
        int u, v;
        scanf("%d %d %lld", &u, &v, &weight[i]);
        g.add_edge(u, v, weight[i]);
    }

    int k;
    scanf("%d", &k);
    vector<int> terminals(k);
    for (int& t : terminals) scanf("%d", &t);

    vector<int> edges = g.tree(terminals);
    long long total = 0;
    for (int id : edges) total += weight[id];
    assert(total == g.solve(terminals));

    printf("%lld %d\n", total, (int)edges.size());
    for (int id : edges) printf("%d ", id);
    printf("\n");
    return 0;
}
