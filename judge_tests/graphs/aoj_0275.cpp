// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0275
// competitive-verifier: TLE 8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/reachability.cpp"
#undef main

vector<long long> shortest(const vector<vector<pair<int, int>>>& adj, int src){
    vector<long long> dist(adj.size(), LLONG_MAX / 4);
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
    dist[src] = 0, pq.push({0, src});
    while (!pq.empty()){
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        for (auto [v, w] : adj[u]){
            if (d + w < dist[v]) dist[v] = d + w, pq.push({dist[v], v});
        }
    }

    return dist;
}

int main(){
    int n, m;
    if (scanf("%d %d", &n, &m) != 2) return 0;
    vector<vector<pair<int, int>>> adj(n);
    vector<array<int, 3>> roads(m);
    for (auto& [u, v, w] : roads){
        if (scanf("%d %d %d", &u, &v, &w) != 3) return 0;
        u--, v--;
        adj[u].push_back({v, w}), adj[v].push_back({u, w});
    }

    int a, b, q;
    if (scanf("%d %d %d", &a, &b, &q) != 3) return 0;
    vector<long long> from_a = shortest(adj, a - 1), to_b = shortest(adj, b - 1);
    long long total = from_a[b - 1];

    /// Since c != d, the ticket allows c -> d iff d is reachable over edges lying on shortest a -> b paths
    Reachability g(n);
    for (auto [u, v, w] : roads){
        for (int it = 0; it < 2; it++, swap(u, v)){
            if (from_a[u] + w + to_b[v] == total) g.add_edge(u, v);
        }
    }

    vector<pair<int, int>> queries(q);
    for (auto& [c, d] : queries){
        if (scanf("%d %d", &c, &d) != 2) return 0;
        c--, d--;
    }
    vector<char> res = g.reaches(queries);
    for (char r : res) puts(r ? "Yes" : "No");

    return 0;
}
