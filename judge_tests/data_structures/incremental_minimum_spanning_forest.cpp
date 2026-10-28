// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/incremental_minimum_spanning_forest
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/link_cut_tree.cpp"
#undef main

int main(){
    int n, m;
    scanf("%d %d", &n, &m);

    /// Edges become nodes n..n + m - 1 holding weight * 2^20 + index, distinct weights make the path max name one edge
    auto heaviest = [](long long a, long long b){ return max(a, b); };
    LinkCutTree<long long, decltype(heaviest)> lct(n + m, -1, heaviest);
    vector<int> eu(m), ev(m);
    for (int i = 0; i < m; i++){
        long long w;
        scanf("%d %d %lld", &eu[i], &ev[i], &w);
        int u = eu[i], v = ev[i], e = n + i;

        int res = -1;
        if (lct.connected(u, v)){
            long long top = lct.query(u, v);
            if (top < (w << 20)){
                res = i;
            }
            else{
                res = top & ((1 << 20) - 1);
                lct.cut(n + res, eu[res]);
                lct.cut(n + res, ev[res]);
            }
        }

        if (res != i){
            lct.set(e, w << 20 | i);
            lct.link(e, u);
            lct.link(v, e);
        }
        printf("%d%c", res, i + 1 == m ? '\n' : ' ');
    }

    return 0;
}
