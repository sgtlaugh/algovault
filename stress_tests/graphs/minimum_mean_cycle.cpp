#include "../common.h"

#define main library_main
#include "../../code_library/graphs/minimum_mean_cycle.cpp"
#undef main

using Edges = vector<array<long long, 3>>;

/// Every simple cycle once, rooted at its smallest vertex, returns the minimum (sum, length) or {0, 0} if acyclic
pair<__int128, long long> brute_min_mean(int n, const Edges& edges){
    vector<vector<int>> out(n);
    for (int i = 0; i < (int)edges.size(); i++) out[edges[i][0]].push_back(i);

    pair<__int128, long long> best = {0, 0};
    vector<bool> on_path(n, false);
    function<void(int, int, __int128, long long)> dfs = [&](int s, int u, __int128 sum, long long len){
        for (int e : out[u]){
            int v = edges[e][1];
            __int128 total = sum + edges[e][2];
            if (v == s){
                if (best.second == 0 || total * best.second < best.first * (len + 1)) best = {total, len + 1};
            }
            else if (v > s && !on_path[v]){
                on_path[v] = true;
                dfs(s, v, total, len + 1);
                on_path[v] = false;
            }
        }
    };

    for (int s = 0; s < n; s++){
        on_path[s] = true;
        dfs(s, s, 0, 0);
        on_path[s] = false;
    }
    return best;
}

/// The reported cycle must be a simple closed walk over real edges whose mean is exactly num / den
void check_cycle(int n, const Edges& edges, const MinimumMeanCycle& g){
    int len = g.cycle.size();
    assert(len >= 1 && len <= n);

    __int128 sum = 0;
    set<int> vertices;
    for (int i = 0; i < len; i++){
        int e = g.cycle[i], next = g.cycle[(i + 1) % len];
        assert(e >= 0 && e < (int)edges.size());
        assert(edges[e][1] == edges[next][0]);
        sum += edges[e][2];
        vertices.insert(edges[e][0]);
    }
    assert((int)vertices.size() == len);
    assert(sum * g.den == (__int128)g.num * len);
}

/// Bellman Ford under weights w * den - num: no negative cycle means no cycle has mean below num / den
bool no_cycle_below(int n, const Edges& edges, long long num, long long den){
    vector<__int128> h(n, 0);
    for (int round = 0; round <= n; round++){
        bool changed = false;
        for (auto [u, v, w] : edges){
            __int128 c = (__int128)w * den - num;
            if (h[u] + c < h[v]) h[v] = h[u] + c, changed = true;
        }
        if (!changed) return true;
    }
    return false;
}

/// Kahn's algorithm removes every vertex iff the graph has no cycle
bool is_acyclic(int n, const Edges& edges){
    vector<int> indeg(n, 0), order;
    vector<vector<int>> out(n);
    for (auto [u, v, w] : edges) out[u].push_back(v), indeg[v]++;

    for (int v = 0; v < n; v++){
        if (indeg[v] == 0) order.push_back(v);
    }
    for (int i = 0; i < (int)order.size(); i++){
        for (int v : out[order[i]]){
            if (--indeg[v] == 0) order.push_back(v);
        }
    }
    return (int)order.size() == n;
}

MinimumMeanCycle build(int n, const Edges& edges){
    MinimumMeanCycle g(n);
    for (auto [u, v, w] : edges) g.add_edge(u, v, w);
    return g;
}

Edges random_edges(int n, int m, long long lo, long long hi){
    Edges edges;
    for (int i = 0; i < m; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), stress::rand_int(lo, hi)});
    return edges;
}

void compare_brute(int n, const Edges& edges, MinimumMeanCycle& g){
    bool ok = g.solve();
    auto [sum, len] = brute_min_mean(n, edges);
    assert(ok == (len != 0));
    if (!ok){
        assert(g.cycle.empty());
        return;
    }

    assert(g.den > 0 && gcd(g.num, g.den) == 1);
    assert(sum * g.den == (__int128)g.num * len);
    check_cycle(n, edges, g);
}

void compare_brute(int n, const Edges& edges){
    MinimumMeanCycle g = build(n, edges);
    compare_brute(n, edges, g);
}

int main(){
    assert(!MinimumMeanCycle(0).solve());

    for (long long it = 0; it < stress::scaled(6000); it++){
        int n = stress::rand_int(1, 7), m = stress::rand_int(0, 2 * n + 2);
        long long range = it % 3 == 0 ? 3 : 1000;
        compare_brute(n, random_edges(n, m, -range, range));
    }

    /// One instance re-solved after every add_edge, so state left by an earlier solve must not leak
    for (long long it = 0; it < stress::scaled(1000); it++){
        int n = stress::rand_int(1, 6), m = stress::rand_int(1, 2 * n + 2);
        Edges edges = random_edges(n, m, -50, 50), prefix;
        MinimumMeanCycle g(n);
        for (auto [u, v, w] : edges){
            g.add_edge(u, v, w), prefix.push_back({u, v, w});
            compare_brute(n, prefix, g);
        }
    }

    /// Weights at the documented bound n * |w| <= 4e18, the brute force sums in __int128
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(1, 4), m = stress::rand_int(1, 7);
        long long bound = 4000000000000000000LL / n;
        Edges edges = random_edges(n, m, -bound, bound);
        for (auto& e : edges){
            int pick = stress::rand_int(0, 2);
            if (pick < 2) e[2] = pick == 0 ? bound : -bound;
        }
        compare_brute(n, edges);
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(2, 150), m = stress::rand_int(n, 6 * n);
        long long range = stress::rand_int(1, 1000000);
        Edges edges = random_edges(n, m, -range, range);
        MinimumMeanCycle g = build(n, edges);
        if (!g.solve()){
            assert(is_acyclic(n, edges));
            continue;
        }
        check_cycle(n, edges, g);
        assert(no_cycle_below(n, edges, g.num, g.den));
    }

    /// One long cycle plus heavy chords: the optimum uses all n edges, so the walk table is fully exercised
    int n = 1500;
    Edges edges;
    for (int v = 0; v < n; v++) edges.push_back({v, (v + 1) % n, v == 0 ? 1 : 0});
    for (int i = 0; i < 3 * n; i++) edges.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1), 1000000});
    MinimumMeanCycle g = build(n, edges);
    assert(g.solve() && g.num == 1 && g.den == n);
    check_cycle(n, edges, g);
    assert(no_cycle_below(n, edges, g.num, g.den));

    return 0;
}
