// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2251
// competitive-verifier: TLE 8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/minimum_path_cover.cpp"
#undef main

int main(){
    const long long INF = LLONG_MAX / 4;
    int n, m, l;
    while (scanf("%d %d %d", &n, &m, &l) == 3 && (n || m || l)){
        vector<vector<long long>> dist(n, vector<long long>(n, INF));
        for (int i = 0; i < n; i++) dist[i][i] = 0;
        for (int i = 0; i < m; i++){
            int u, v;
            long long d;
            if (scanf("%d %d %lld", &u, &v, &d) != 3) return 0;
            dist[u][v] = dist[v][u] = min(dist[u][v], d);
        }
        for (int k = 0; k < n; k++){
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++) dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }

        vector<pair<long long, int>> requests(l);
        for (auto& [t, p] : requests){
            if (scanf("%d %lld", &p, &t) != 2) return 0;
        }
        sort(requests.begin(), requests.end());

        /// Edges only go forward in time order, so ties at distance 0 cannot close a cycle
        DAGPathCover g(l);
        for (int i = 0; i < l; i++){
            for (int j = i + 1; j < l; j++){
                auto [ti, pi] = requests[i];
                auto [tj, pj] = requests[j];
                if (ti + dist[pi][pj] <= tj) g.add_edge(i, j);
            }
        }

        int cover = g.min_disjoint_path_cover();
        assert(g.min_path_cover() == cover && (int)g.max_antichain().size() == cover);
        printf("%d\n", cover);
    }
    return 0;
}
