// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/eulerian_trail_undirected
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/euler_path.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;

    string out;
    while (t--){
        int n, m;
        if (scanf("%d %d", &n, &m) != 2) return 0;

        EulerUndirected g(n);
        for (int i = 0; i < m; i++){
            int u, v;
            if (scanf("%d %d", &u, &v) != 2) return 0;
            g.add_edge(u, v);
        }

        auto path = g.solve();
        if (path.empty()){
            out += "No\n";
            continue;
        }

        out += "Yes\n";
        for (int i = 0; i <= m; i++) out += to_string(path[i]) + (i < m ? " " : "\n");
        for (int i = 0; i < m; i++) out += to_string(g.walk_edges[i]) + (i + 1 < m ? " " : "");
        out += "\n";
    }
    fputs(out.c_str(), stdout);
    return 0;
}
