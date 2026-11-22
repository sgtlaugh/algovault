// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2171
// competitive-verifier: TLE 8
// competitive-verifier: ERROR 1e-8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/gauss_jordan.cpp"
#undef main

int main(){
    int n, s, t;
    while (scanf("%d %d %d", &n, &s, &t) == 3 && n){
        s--, t--;
        vector<int> sign(n);
        vector<vector<int>> a(n, vector<int>(n));
        for (int& x : sign) if (scanf("%d", &x) != 1) return 0;
        for (auto& row : a){
            for (int& x : row) if (scanf("%d", &x) != 1) return 0;
        }

        const int INF = INT_MAX / 2;
        vector<int> dist(n, INF);
        vector<char> done(n, 0);
        dist[t] = 0;
        for (int it = 0; it < n; it++){
            int u = -1;
            for (int v = 0; v < n; v++){
                if (!done[v] && (u == -1 || dist[v] < dist[u])) u = v;
            }
            if (dist[u] == INF) break;
            done[u] = 1;
            for (int v = 0; v < n; v++){
                if (a[u][v]) dist[v] = min(dist[v], dist[u] + a[u][v]);
            }
        }

        if (dist[s] == INF){
            puts("impossible");
            continue;
        }

        /// E[u] = mean over u's choices v of a[u][v] + E[v], E[t] = 0, over the vertices that reach t
        vector<vector<double>> equations(n, vector<double>(n + 1, 0.0));
        for (int u = 0; u < n; u++){
            equations[u][u] = 1;
            if (u == t || dist[u] == INF) continue;

            vector<int> choices;
            for (int v = 0; v < n; v++){
                if (a[u][v] && (!sign[u] || dist[u] == a[u][v] + dist[v])) choices.push_back(v);
            }
            for (int v : choices){
                equations[u][v] -= 1.0 / choices.size();
                equations[u][n] += (double)a[u][v] / choices.size();
            }
        }

        GaussSolution<double> res = gauss(equations);
        printf("%.10f\n", res.x[s]);
    }

    return 0;
}
