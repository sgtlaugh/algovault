// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/persistent_unionfind
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/disjoint_set.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 0;

    /// Version i + 1 is the graph after query i, version 0 is the empty graph
    vector<int> t(q), u(q), v(q), answer(q, -1);
    vector<vector<int>> children(q + 1), asks(q + 1);
    for (int i = 0; i < q; i++){
        int k;
        if (scanf("%d %d %d %d", &t[i], &k, &u[i], &v[i]) != 4) return 0;
        if (t[i] == 0) children[k + 1].push_back(i + 1);
        else asks[k + 1].push_back(i);
    }

    RollbackDSU dsu(n);
    vector<int> entry_time(q + 1);
    vector<pair<int, bool>> stack = {{0, false}};
    while (!stack.empty()){
        auto [node, leaving] = stack.back();
        stack.pop_back();
        if (leaving){
            dsu.rollback(entry_time[node]);
            continue;
        }

        entry_time[node] = dsu.time();
        if (node > 0) dsu.connect(u[node - 1], v[node - 1]);
        for (int i: asks[node]) answer[i] = dsu.is_connected(u[i], v[i]);

        stack.push_back({node, true});
        for (int child: children[node]) stack.push_back({child, false});
    }

    for (int i = 0; i < q; i++){
        if (t[i] == 1) printf("%d\n", answer[i]);
    }
    return 0;
}
