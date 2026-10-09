#include "../common.h"

#define main library_main
#include "../../code_library/graphs/k_shortest_walks.cpp"
#undef main

using Edges = vector<array<long long, 3>>;

/// Enumerates walk prefixes best-first, needs weights >= 1: a zero weight cycle before t would be popped forever
/// Parallel cheap loops make the prefix count explode, so it gives up (complete = false) after 3e6 pushes
vector<long long> enumerate_walks(int n, const Edges& edges, int s, int t, int k, bool& complete){
    vector<char> reaches_t(n, 0);
    reaches_t[t] = 1;
    for (int round = 0; round < n; round++){
        for (auto [u, v, w] : edges){
            if (reaches_t[v]) reaches_t[u] = 1;
        }
    }

    vector<long long> result;
    complete = true;
    if (!reaches_t[s] || k <= 0) return result;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> prefixes;
    long long pushed = 1;
    prefixes.push({0, s});
    while ((int)result.size() < k && !prefixes.empty()){
        auto [len, u] = prefixes.top();
        prefixes.pop();
        if (u == t) result.push_back(len);
        for (auto [x, v, w] : edges){
            if (x == u && reaches_t[v]) prefixes.push({len + w, (int)v}), pushed++;
        }
        if (pushed >= 3000000){
            complete = false;
            break;
        }
    }
    return result;
}

/// Dijkstra that settles every vertex up to k times, the j-th settle of t is the j-th shortest walk
vector<long long> k_settle_dijkstra(int n, const Edges& edges, int s, int t, int k){
    vector<vector<pair<int, long long>>> adj(n);
    for (auto [u, v, w] : edges) adj[u].push_back({(int)v, w});

    vector<int> settled(n, 0);
    vector<long long> result;
    priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> heap;
    heap.push({0, s});
    while ((int)result.size() < k && !heap.empty()){
        auto [len, u] = heap.top();
        heap.pop();
        if (settled[u] >= k) continue;
        settled[u]++;
        if (u == t) result.push_back(len);
        for (auto [v, w] : adj[u]) heap.push({len + w, v});
    }
    return result;
}

Edges random_edges(int n, int m, long long min_w, long long max_w){
    Edges edges;
    for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(min_w, max_w)});
    return edges;
}

/// Path v -> v + 1 with a loop of weight n - v on v: a walk is the path plus one partition of its extra length into parts <= n
Edges loop_ladder(int n){
    Edges edges;
    for (int v = 0; v < n; v++){
        if (v + 1 < n) edges.push_back({v, v + 1, 1});
        edges.push_back({v, v, n - v});
    }
    return edges;
}

vector<long long> loop_ladder_walks(int n, int k){
    vector<long long> result;
    for (int extra = 0; (int)result.size() < k; extra++){
        vector<long long> partitions(extra + 1, 0);
        partitions[0] = 1;
        for (int part = 1; part <= min(n, extra); part++){
            for (int j = part; j <= extra; j++) partitions[j] += partitions[j - part];
        }
        for (long long c = 0; c < partitions[extra] && (int)result.size() < k; c++) result.push_back(n - 1 + extra);
    }
    return result;
}

vector<long long> run_library(int n, const Edges& edges, int s, int t, int k){
    KShortestWalks g(n);
    for (auto [u, v, w] : edges) g.add_edge(u, v, w);
    return g.solve(s, t, k);
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 5), m = stress::rand_int(0, 9);
        Edges edges = random_edges(n, m, 1, it % 3 == 0 ? 2 : 5);
        int s = stress::rand_int(0, n - 1), t = stress::rand_int(0, n - 1), k = stress::rand_int(0, 25);
        bool complete;
        vector<long long> walks = enumerate_walks(n, edges, s, t, k, complete);
        assert(run_library(n, edges, s, t, k) == (complete ? walks : k_settle_dijkstra(n, edges, s, t, k)));

        Edges with_zeros = random_edges(n, m, 0, 2);
        assert(run_library(n, with_zeros, s, t, k) == k_settle_dijkstra(n, with_zeros, s, t, k));
    }

    for (long long it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(2, 60), m = stress::rand_int(0, 4 * n);
        Edges edges = random_edges(n, m, 0, it % 2 ? 1000000000 : 3);
        int s = stress::rand_int(0, n - 1), t = stress::rand_int(0, n - 1), k = stress::rand_int(1, 300);

        /// several queries on one graph check that solve leaves no state behind
        KShortestWalks g(n);
        for (auto [u, v, w] : edges) g.add_edge(u, v, w);
        assert(g.solve(s, t, k) == k_settle_dijkstra(n, edges, s, t, k));
        assert(g.solve(t, s, k) == k_settle_dijkstra(n, edges, t, s, k));
    }

    /// Library Checker scale: 3e5 edges with weights up to 1e9
    int n = 100000, m = 300000;
    Edges edges = random_edges(n, m, 0, 1000000000);
    assert(run_library(n, edges, 0, n - 1, 5) == k_settle_dijkstra(n, edges, 0, n - 1, 5));

    /// A path with a weight 1 loop on each of its n vertices: one walk of length n - 1, n of length n, then C(n + 1, 2) of length n + 1
    n = 200000;
    edges.clear();
    for (int v = 0; v < n; v++){
        if (v + 1 < n) edges.push_back({v, v + 1, 1});
        edges.push_back({v, v, 1});
    }
    vector<long long> walks = run_library(n, edges, 0, n - 1, 300000);
    assert((int)walks.size() == 300000 && walks[0] == n - 1);
    assert(count(walks.begin(), walks.end(), n) == n && count(walks.begin(), walks.end(), n + 1) == 300000 - 1 - n);

    /// Loops of growing weight push every new sidetrack to the bottom of its heap, only rank balancing keeps the melds O(log m)
    n = 3000;
    KShortestWalks ladder(n);
    for (auto [u, v, w] : loop_ladder(n)) ladder.add_edge(u, v, w);
    assert(ladder.solve(0, n - 1, 20000) == loop_ladder_walks(n, 20000));
    assert(ladder.pool.size() <= ladder.edges.size() + 2 * n * (__lg(ladder.edges.size()) + 2));

    n = 200000;
    assert(run_library(n, loop_ladder(n), 0, n - 1, 300000) == loop_ladder_walks(n, 300000));
    return 0;
}
