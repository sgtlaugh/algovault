// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/1605
// competitive-verifier: TLE 8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/weighted_matroid_intersection.cpp"
#undef main

int main(){
    int n, m, k;
    while (scanf("%d %d %d", &n, &m, &k) == 3 && (n || m || k)){
        GraphicMatroid forest(n);
        vector<int> company(m);
        vector<long long> cost(m);
        for (int i = 0; i < m; i++){
            int u, v;
            char c;
            if (scanf("%d %d %lld %c", &u, &v, &cost[i], &c) != 4) return 0;
            forest.add_edge(u - 1, v - 1);
            company[i] = c == 'A' ? 0 : 1;
        }

        /// Caps k and n - 1 - k make every n - 1 edge common independent set a spanning tree with exactly k A edges
        PartitionMatroid split(company, {k, n - 1 - k});

        /// Every weight gets B more than any cost sum, so the largest sets win first and among them the cheapest
        long long big = 2 * accumulate(cost.begin(), cost.end(), 0LL) + 1;
        vector<long long> weight(m);
        for (int i = 0; i < m; i++) weight[i] = big - cost[i];

        vector<int> chosen = weighted_matroid_intersection(m, weight, forest, split);
        long long total = 0;
        for (int x : chosen) total += cost[x];
        printf("%lld\n", (int)chosen.size() == n - 1 ? total : -1);
    }

    return 0;
}
