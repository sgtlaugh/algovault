#include "../common.h"

#define main library_main
#include "../../code_library/graphs/find_cycle.cpp"
#undef main

/// Kahn's algorithm when directed, union-find when undirected: a cycle exists iff they fail
bool brute_has_cycle(int n, const vector<pair<int, int>>& es, bool directed){
    if (directed){
        vector<int> indegree(n, 0), queue;
        vector<vector<int>> out(n);
        for (auto [u, v] : es) out[u].push_back(v), indegree[v]++;
        for (int v = 0; v < n; v++) if (!indegree[v]) queue.push_back(v);
        for (int i = 0; i < (int)queue.size(); i++){
            for (int v : out[queue[i]]) if (!--indegree[v]) queue.push_back(v);
        }

        return (int)queue.size() < n;
    }

    vector<int> parent(n);
    iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x){
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (auto [u, v] : es){
        if (find(u) == find(v)) return true;
        parent[find(u)] = find(v);
    }

    return false;
}

void check_cycle(int n, const vector<pair<int, int>>& es, bool directed, const FindCycle::Cycle& c){
    int k = c.vertices.size();
    assert(k >= 1 && (int)c.edges.size() == k);
    assert((int)set<int>(c.vertices.begin(), c.vertices.end()).size() == k);
    assert((int)set<int>(c.edges.begin(), c.edges.end()).size() == k);

    for (int i = 0; i < k; i++){
        int a = c.vertices[i], b = c.vertices[(i + 1) % k], id = c.edges[i];
        assert(0 <= a && a < n && 0 <= id && id < (int)es.size());
        auto [u, v] = es[id];
        assert((u == a && v == b) || (!directed && u == b && v == a));
    }
}

void check(int n, const vector<pair<int, int>>& es, bool directed){
    FindCycle g(n, directed ? FindCycle::DIRECTED : FindCycle::UNDIRECTED);
    for (auto [u, v] : es) g.add_edge(u, v);
    FindCycle::Cycle c = g.solve();

    assert(c.exists == brute_has_cycle(n, es, directed));
    if (c.exists) check_cycle(n, es, directed, c);
    else assert(c.vertices.empty() && c.edges.empty());
}

vector<pair<int, int>> random_edges(int n, int m, int shape){
    vector<pair<int, int>> es;
    /// Shape 1 keeps edges pointing up a random order, a DAG or a forest-heavy graph, so "no cycle" is common
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), stress::rng());
    for (int i = 0; i < m; i++){
        int a = stress::rand_int(0, n - 1), b = stress::rand_int(0, n - 1);
        if (shape == 1 && a == b) continue;
        if (shape == 1 && a > b) swap(a, b);
        es.push_back(shape == 1 ? make_pair(order[a], order[b]) : make_pair(a, b));
    }

    return es;
}

int main(){
    /// Every graph with up to 3 vertices and 3 edges, self loops and parallel edges included, both kinds
    for (int n = 1; n <= 3; n++){
        int pairs = n * n;
        for (int m = 0; m <= 3; m++){
            int total = 1;
            for (int i = 0; i < m; i++) total *= pairs;
            for (int code = 0; code < total; code++){
                vector<pair<int, int>> es;
                for (int i = 0, x = code; i < m; i++, x /= pairs) es.push_back({x % pairs / n, x % pairs % n});
                check(n, es, true);
                check(n, es, false);
            }
        }
    }

    for (int it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(1, 40), m = stress::rand_int(0, 60);
        vector<pair<int, int>> es = random_edges(n, m, it % 2);
        check(n, es, true);
        check(n, es, false);
    }

    /// Forests: an undirected spanning tree plus isolated vertices has no cycle until one more edge joins a tree
    for (int it = 0; it < stress::scaled(300); it++){
        int n = stress::rand_int(2, 200);
        vector<pair<int, int>> es;
        for (int v = 1; v < n; v++) if (stress::rand_int(0, 4)) es.push_back({stress::rand_int(0, v - 1), v});
        check(n, es, false);
        es.push_back({stress::rand_int(0, n - 1), stress::rand_int(0, n - 1)});
        check(n, es, false);
    }

    /// A path of 5e5 vertices closed into one cycle, both kinds: the DFS stack is as deep as the graph
    int n = 500000;
    vector<pair<int, int>> ring;
    for (int v = 0; v + 1 < n; v++) ring.push_back({v, v + 1});
    check(n, ring, true);
    check(n, ring, false);
    ring.push_back({n - 1, 0});
    check(n, ring, true);
    check(n, ring, false);

    stress::expect_assert("graphs/find_cycle.cpp", []{ FindCycle bad(2, FindCycle::DIRECTED); bad.add_edge(0, 2); });

    return 0;
}
