#include "../common.h"

#define main library_main
#include "../../code_library/graphs/gomory_hu.cpp"
#undef main

struct Link{ int u, v; long long cap; };

/// Min cut by enumerating every side containing s but not t
long long subset_min_cut(int n, int s, int t, const vector<Link>& links){
    long long best = LLONG_MAX;
    for (int mask = 0; mask < (1 << n); mask++){
        if (!(mask >> s & 1) || (mask >> t & 1)) continue;
        long long cut = 0;
        for (auto& e : links) if ((mask >> e.u & 1) != (mask >> e.v & 1)) cut += e.cap;
        best = min(best, cut);
    }
    return best;
}

/// Edmonds-Karp on a capacity matrix, a max flow independent of the library's Dinic
long long matrix_max_flow(int n, int s, int t, const vector<Link>& links){
    vector<vector<long long>> cap(n, vector<long long>(n, 0));
    for (auto& e : links){
        if (e.u != e.v) cap[e.u][e.v] += e.cap, cap[e.v][e.u] += e.cap;
    }

    long long flow = 0;
    while (true){
        vector<int> from(n, -1);
        vector<int> queue = {s};
        from[s] = s;
        for (int i = 0; i < (int)queue.size() && from[t] == -1; i++){
            int u = queue[i];
            for (int v = 0; v < n; v++){
                if (from[v] == -1 && cap[u][v] > 0) from[v] = u, queue.push_back(v);
            }
        }
        if (from[t] == -1) return flow;

        long long push = LLONG_MAX;
        for (int v = t; v != s; v = from[v]) push = min(push, cap[from[v]][v]);
        for (int v = t; v != s; v = from[v]) cap[from[v]][v] -= push, cap[v][from[v]] += push;
        flow += push;
    }
}

vector<Link> random_links(int n, int m, long long max_cap){
    vector<Link> links;
    for (int i = 0; i < m; i++){
        links.push_back({(int)stress::rand_int(0, n - 1), (int)stress::rand_int(0, n - 1), stress::rand_int(0, max_cap)});
    }
    return links;
}

GomoryHu build_tree(int n, const vector<Link>& links){
    GomoryHu g(n);
    for (auto& e : links) g.add_edge(e.u, e.v, e.cap);
    g.build();
    return g;
}

/// Every pair against a reference min cut, each tree edge carries the cut between its endpoints,
/// and the subtree below each tree edge is a side of the graph whose cut equals the edge weight
template<typename Reference>
void check_all_pairs(int n, const vector<Link>& links, Reference reference){
    GomoryHu g = build_tree(n, links);
    vector<vector<long long>> cut(n, vector<long long>(n));
    for (int s = 0; s < n; s++){
        for (int t = s + 1; t < n; t++) cut[s][t] = cut[t][s] = reference(n, s, t, links);
    }

    vector<vector<int>> below(n, vector<int>(n, 0));
    for (int x = 0; x < n; x++){
        int u = x, steps = 0;
        for (; u != 0; u = g.tree_parent[u], steps++){
            assert(0 <= g.tree_parent[u] && g.tree_parent[u] < n && steps < n);
            below[u][x] = 1;
        }
    }

    for (int v = 1; v < n; v++){
        assert(g.tree_weight[v] == cut[v][g.tree_parent[v]]);
        long long side_cut = 0;
        for (auto& e : links) if (below[v][e.u] != below[v][e.v]) side_cut += e.cap;
        assert(side_cut == g.tree_weight[v]);
    }
    for (int s = 0; s < n; s++){
        for (int t = 0; t < n; t++) if (s != t) assert(g.min_cut(s, t) == cut[s][t]);
    }
}

int main(){
    for (long long it = 0; it < stress::scaled(2000); it++){
        int n = stress::rand_int(1, 8), m = stress::rand_int(0, 16);
        long long max_cap = it % 4 == 0 ? 1000000000000000LL : it % 4 == 1 ? LLONG_MAX / 16 : 10;
        check_all_pairs(n, random_links(n, m, max_cap), subset_min_cut);
    }

    for (long long it = 0; it < stress::scaled(60); it++){
        int n = stress::rand_int(9, 40), m = stress::rand_int(n - 1, 4 * n);
        long long max_cap = it % 3 ? 20 : 1000000000LL;
        check_all_pairs(n, random_links(n, m, max_cap), matrix_max_flow);
    }

    /// A path with shuffled labels: the min cut is the lightest link between the endpoints, and Dinic recurses about n deep
    {
        int n = 5000;
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), stress::rng());

        vector<Link> links;
        for (int i = 0; i + 1 < n; i++) links.push_back({order[i], order[i + 1], stress::rand_int(1, 1000000000000000LL)});
        GomoryHu g = build_tree(n, links);

        for (int q = 0; q < 2000; q++){
            int i = stress::rand_int(0, n - 1), j = stress::rand_int(0, n - 1);
            if (i == j) continue;
            long long lightest = LLONG_MAX;
            for (int k = min(i, j); k < max(i, j); k++) lightest = min(lightest, links[k].cap);
            assert(g.min_cut(order[i], order[j]) == lightest);
        }
    }

    /// A single edge at the capacity limit
    GomoryHu widest = build_tree(2, {{0, 1, LLONG_MAX}});
    assert(widest.min_cut(0, 1) == LLONG_MAX && widest.tree_weight[1] == LLONG_MAX);
    return 0;
}
