#include "../common.h"

#define main library_main
#include "../../code_library/graphs/global_min_cut.cpp"
#undef main

long long cut_weight(const vector<vector<long long>>& w, const vector<bool>& in_side){
    long long total = 0;
    for (int u = 0; u < (int)w.size(); u++){
        for (int v = u + 1; v < (int)w.size(); v++){
            if (in_side[u] != in_side[v]) total += w[u][v];
        }
    }
    return total;
}

/// Every cut, node 0 fixed outside the side so each of the 2^(n-1) - 1 cuts is seen once
long long brute_min_cut(const vector<vector<long long>>& w){
    int n = w.size();
    long long best = LLONG_MAX;
    for (int mask = 1; mask < (1 << (n - 1)); mask++){
        vector<bool> in_side(n, false);
        for (int i = 1; i < n; i++) in_side[i] = mask >> (i - 1) & 1;
        best = min(best, cut_weight(w, in_side));
    }
    return best;
}

/// Global min cut = min over t of the 0-t min cut, each found with Edmonds-Karp on the matrix
long long flow_min_cut(const vector<vector<long long>>& w){
    int n = w.size();
    long long best = LLONG_MAX;
    for (int t = 1; t < n; t++){
        auto cap = w;
        long long flow = 0;
        while (true){
            vector<int> par(n, -1);
            queue<int> q;
            par[0] = 0, q.push(0);
            while (!q.empty() && par[t] == -1){
                int u = q.front();
                q.pop();
                for (int v = 0; v < n; v++){
                    if (par[v] == -1 && cap[u][v] > 0) par[v] = u, q.push(v);
                }
            }
            if (par[t] == -1) break;

            long long push = LLONG_MAX;
            for (int v = t; v != 0; v = par[v]) push = min(push, cap[par[v]][v]);
            for (int v = t; v != 0; v = par[v]) cap[par[v]][v] -= push, cap[v][par[v]] += push;
            flow += push;
        }
        best = min(best, flow);
    }
    return best;
}

/// The side must be a non-empty proper subset of distinct nodes whose cut weighs exactly the returned weight
void check_side(const vector<vector<long long>>& w, long long weight, const vector<int>& side){
    int n = w.size();
    assert(!side.empty() && (int)side.size() < n);
    vector<bool> in_side(n, false);
    for (int x: side){
        assert(0 <= x && x < n && !in_side[x]);
        in_side[x] = true;
    }
    assert(cut_weight(w, in_side) == weight);
}

/// Random weights, sparse to dense, small weights for ties and huge ones summing close to LLONG_MAX
pair<GlobalMinCut, vector<vector<long long>>> random_graph(int n, int edges, long long max_w){
    GlobalMinCut g(n);
    vector<vector<long long>> w(n, vector<long long>(n, 0));
    for (int e = 0; e < edges; e++){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        long long weight = stress::rand_int(0, max_w);
        g.add_edge(u, v, weight);
        if (u != v) w[u][v] += weight, w[v][u] += weight;
    }
    return {g, w};
}

int main(){
    for (long long it = 0; it < stress::scaled(3000); it++){
        int n = stress::rand_int(2, 12), edges = stress::rand_int(0, n * n);
        long long max_w = it % 3 == 0 ? 3 : it % 3 == 1 ? 1000000 : (LLONG_MAX - 1) / max(1, edges);

        auto [g, w] = random_graph(n, edges, max_w);
        auto [weight, side] = g.min_cut();
        assert(weight == brute_min_cut(w));
        check_side(w, weight, side);

        /// min_cut works on a copy, so a query after adding more edges sees the whole graph
        if (it % 3 == 2) continue;
        int u = stress::rand_int(0, n - 1), v = (u + stress::rand_int(1, n - 1)) % n;
        long long extra = stress::rand_int(0, max_w);
        g.add_edge(u, v, extra), w[u][v] += extra, w[v][u] += extra;
        tie(weight, side) = g.min_cut();
        assert(weight == brute_min_cut(w));
        check_side(w, weight, side);
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        int n = stress::rand_int(13, 40), edges = stress::rand_int(0, 2 * n * n);
        auto [g, w] = random_graph(n, edges, it % 2 ? 2 : 1000000000);

        auto [weight, side] = g.min_cut();
        assert(weight == flow_min_cut(w));
        check_side(w, weight, side);
    }

    /// Total weight of exactly LLONG_MAX: at n = 2 the only cut weighs LLONG_MAX, the starting best
    for (long long split: {0LL, 5LL, LLONG_MAX / 2}){
        GlobalMinCut g(2);
        vector<vector<long long>> w = {{0, LLONG_MAX}, {LLONG_MAX, 0}};
        g.add_edge(0, 1, LLONG_MAX - split), g.add_edge(1, 0, split);
        auto [weight, side] = g.min_cut();
        assert(weight == LLONG_MAX);
        check_side(w, weight, side);
    }

    /// Complete graph with unit weights at the O(n^3) worst case: every single node is a cut of n - 1
    {
        int n = 400;
        GlobalMinCut g(n);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) g.add_edge(u, v, 1);
        }
        auto [weight, side] = g.min_cut();
        assert(weight == n - 1 && (side.size() == 1 || (int)side.size() == n - 1));
    }

    return 0;
}
