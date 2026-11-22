#include "../common.h"

#define main library_main
#include "../../code_library/graphs/minimum_weight_cycle.cpp"
#undef main

struct RawEdge{
    int u, v;
    long long w;
};

/// Every simple cycle once per starting vertex, the smallest vertex on it, each edge id used at most once
long long brute_cycles(int n, const vector<RawEdge>& es, bool directed, int s, int cur, vector<char>& seen, vector<char>& used, long long sum){
    long long best = -1;
    for (int id = 0; id < (int)es.size(); id++){
        if (used[id]) continue;
        int to;
        if (es[id].u == cur) to = es[id].v;
        else if (!directed && es[id].v == cur) to = es[id].u;
        else continue;

        long long cand = -1;
        if (to == s) cand = sum + es[id].w;
        else if (to > s && !seen[to]){
            seen[to] = used[id] = 1;
            cand = brute_cycles(n, es, directed, s, to, seen, used, sum + es[id].w);
            seen[to] = used[id] = 0;
        }
        if (cand != -1 && (best == -1 || cand < best)) best = cand;
    }

    return best;
}

long long brute(int n, const vector<RawEdge>& es, bool directed){
    long long best = -1;
    for (int s = 0; s < n; s++){
        vector<char> seen(n, 0), used(es.size(), 0);
        seen[s] = 1;
        long long cand = brute_cycles(n, es, directed, s, s, seen, used, 0);
        if (cand != -1 && (best == -1 || cand < best)) best = cand;
    }

    return best;
}

/// The lightest cycle through edge u -> v is w plus the shortest v -> u path avoiding that edge
long long edge_removal(int n, const vector<RawEdge>& es, bool directed){
    long long best = -1;
    for (int skip = 0; skip < (int)es.size(); skip++){
        vector<long long> dist(n, LLONG_MAX);
        priority_queue<pair<long long, int>, vector<pair<long long, int>>, greater<>> pq;
        dist[es[skip].v] = 0;
        pq.push({0, es[skip].v});
        while (!pq.empty()){
            auto [d, x] = pq.top();
            pq.pop();
            if (d != dist[x]) continue;
            for (int id = 0; id < (int)es.size(); id++){
                if (id == skip) continue;
                int to;
                if (es[id].u == x) to = es[id].v;
                else if (!directed && es[id].v == x) to = es[id].u;
                else continue;
                if (d + es[id].w < dist[to]) dist[to] = d + es[id].w, pq.push({dist[to], to});
            }
        }
        if (dist[es[skip].u] == LLONG_MAX) continue;
        long long cand = dist[es[skip].u] + es[skip].w;
        if (best == -1 || cand < best) best = cand;
    }

    return best;
}

void check_cycle(int n, const vector<RawEdge>& es, bool directed, const WeightCycle& c, long long expected){
    assert(c.exists == (expected != -1));
    if (!c.exists){
        assert(c.weight == 0 && c.vertices.empty() && c.edges.empty());
        return;
    }
    assert(c.weight == expected);

    int k = c.vertices.size();
    assert(k >= 1 && (int)c.edges.size() == k);
    vector<char> seen(n, 0), used(es.size(), 0);
    long long sum = 0;
    for (int i = 0; i < k; i++){
        int a = c.vertices[i], b = c.vertices[(i + 1) % k], id = c.edges[i];
        assert(!seen[a] && !used[id]);
        seen[a] = used[id] = 1;
        bool forward = es[id].u == a && es[id].v == b;
        bool backward = es[id].u == b && es[id].v == a;
        assert(forward || (!directed && backward));
        sum += es[id].w;
    }
    assert(sum == expected);
}

WeightCycle run(int n, const vector<RawEdge>& es, bool directed){
    MinimumWeightCycle g(n, directed ? MinimumWeightCycle::DIRECTED : MinimumWeightCycle::UNDIRECTED);
    for (auto [u, v, w] : es) g.add_edge(u, v, w);

    return g.solve();
}

/// Weight styles: random, small with zeros and ties, all equal to 1, all equal to 0, all equal to a larger value
vector<RawEdge> random_graph(int n, int m, int style, bool loops){
    long long same = style == 2 ? 1 : style == 3 ? 0 : stress::rand_int(2, 1000);
    vector<RawEdge> es;
    while ((int)es.size() < m){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (u == v && !loops && n > 1) continue;

        long long w = style == 0 ? stress::rand_int(0, 1000000000) : style == 1 ? stress::rand_int(0, 3) : same;
        es.push_back({u, v, w});
    }

    return es;
}

int main(){
    for (int directed = 0; directed < 2; directed++){
        assert(!run(0, {}, directed).exists);
        assert(!run(4, {}, directed).exists);
    }

    /// Small graphs against every simple cycle: parallel edges, self loops, zero weights, isolated vertices
    for (int it = 0; it < stress::scaled(30000); it++){
        int n = stress::rand_int(1, 6), m = stress::rand_int(0, 8);
        bool directed = it % 2;
        vector<RawEdge> es = random_graph(n, m, it / 2 % 5, it % 7 != 0);
        check_cycle(n, es, directed, run(n, es, directed), brute(n, es, directed));
    }

    /// Forests have no undirected cycle, DAGs no directed one, even with many parallel paths
    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(1, 40);
        vector<RawEdge> tree, dag;
        for (int v = 1; v < n; v++){
            if (stress::rand_int(0, 4)) tree.push_back({(int)stress::rand_int(0, v - 1), v, stress::rand_int(0, 9)});
        }
        for (int i = 0; i < 3 * n; i++){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            if (u < v) dag.push_back({u, v, stress::rand_int(0, 9)});
        }
        check_cycle(n, tree, false, run(n, tree, false), -1);
        check_cycle(n, dag, true, run(n, dag, true), -1);
    }

    /// Medium graphs against removing each edge in turn
    for (int it = 0; it < stress::scaled(400); it++){
        int n = stress::rand_int(2, 40), m = stress::rand_int(0, 2 * n);
        bool directed = it % 2;
        vector<RawEdge> es = random_graph(n, m, it / 2 % 5, it % 3 == 0);
        check_cycle(n, es, directed, run(n, es, directed), edge_removal(n, es, directed));
    }

    /// One cycle of 1000 vertices, unweighted and then with its weight near the 1e18 path limit
    int n = 1000;
    vector<RawEdge> ring;
    for (int v = 0; v < n; v++) ring.push_back({v, (v + 1) % n, 1});
    check_cycle(n, ring, false, run(n, ring, false), n);
    check_cycle(n, ring, true, run(n, ring, true), n);
    for (auto& e : ring) e.w = (long long)1e15;
    ring[0].w = 0;
    check_cycle(n, ring, true, run(n, ring, true), (long long)1e15 * (n - 1));

    /// Paths below 1e18 still close cycles up to 2e18: two parallel edges, a triangle, a pair just under the bound
    const long long E17 = 100000000000000000LL, LIMIT = 1000000000000000000LL - 1;
    for (int directed = 0; directed < 2; directed++){
        vector<RawEdge> pair_of = {{0, 1, 6 * E17}, {1, 0, 6 * E17}}, edge_pair = {{0, 1, LIMIT}, {1, 0, LIMIT}};
        vector<RawEdge> triangle = {{0, 1, 4 * E17}, {1, 2, 4 * E17}, {2, 0, 4 * E17}};
        check_cycle(2, pair_of, directed, run(2, pair_of, directed), 12 * E17);
        check_cycle(3, triangle, directed, run(3, triangle, directed), 12 * E17);
        check_cycle(2, edge_pair, directed, run(2, edge_pair, directed), 2 * LIMIT);
    }

    /// A second solve() on the same graph reuses the scratch arrays and must agree with the first
    for (int it = 0; it < stress::scaled(200); it++){
        int n = stress::rand_int(1, 30), m = stress::rand_int(0, 40);
        bool directed = it % 2;
        vector<RawEdge> es = random_graph(n, m, it / 2 % 5, true);
        MinimumWeightCycle g(n, directed ? MinimumWeightCycle::DIRECTED : MinimumWeightCycle::UNDIRECTED);
        for (auto [u, v, w] : es) g.add_edge(u, v, w);
        long long expected = edge_removal(n, es, directed);
        check_cycle(n, es, directed, g.solve(), expected);
        check_cycle(n, es, directed, g.solve(), expected);
    }

    stress::expect_assert("graphs/minimum_weight_cycle.cpp", []{
        MinimumWeightCycle g(2, MinimumWeightCycle::UNDIRECTED);
        g.add_edge(0, 1, -1);
    });

    return 0;
}
