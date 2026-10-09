#include "../common.h"

#define main library_main
#include "../../code_library/graphs/max_clique.cpp"
#undef main

struct Brute{
    int max_clique = 0, max_independent = 0;
    vector<int> maximal;
};

/// Every vertex subset as a bitmask: clique[mask] extends clique[mask without its lowest vertex], and a clique
/// is maximal when no vertex is adjacent to all of it (common = intersection of the neighbourhoods)
Brute brute(int n, const vector<int>& adj){
    Brute result;
    int full = (1 << n) - 1;
    vector<int> common(1 << n), common_complement(1 << n);
    vector<char> clique(1 << n), independent(1 << n);
    common[0] = common_complement[0] = full, clique[0] = independent[0] = 1;

    for (int mask = 1; mask <= full; mask++){
        int v = __builtin_ctz(mask), rest = mask & (mask - 1);
        int complement = full & ~adj[v] & ~(1 << v);
        common[mask] = common[rest] & adj[v];
        common_complement[mask] = common_complement[rest] & complement;
        clique[mask] = clique[rest] && (adj[v] & rest) == rest;
        independent[mask] = independent[rest] && (complement & rest) == rest;
    }

    for (int mask = 0; mask <= full; mask++){
        if (clique[mask]) result.max_clique = max(result.max_clique, __builtin_popcount(mask));
        if (independent[mask]) result.max_independent = max(result.max_independent, __builtin_popcount(mask));
        if (clique[mask] && common[mask] == 0) result.maximal.push_back(mask);
    }
    return result;
}

bool valid(CliqueGraph& g, const vector<int>& vertices, bool independent){
    for (size_t i = 0; i < vertices.size(); i++){
        assert(0 <= vertices[i] && vertices[i] < g.n);
        for (size_t j = 0; j < i; j++){
            if (vertices[i] == vertices[j] || g.adj[vertices[i]][vertices[j]] == independent) return false;
        }
    }
    return true;
}

/// Library answers checked against the brute force: sizes, validity, and the exact set of maximal cliques
void check(CliqueGraph& g, const vector<int>& adj){
    Brute expected = brute(g.n, adj);

    vector<int> clique = g.max_clique(), independent = g.max_independent_set();
    assert((int)clique.size() == expected.max_clique && valid(g, clique, false));
    assert((int)independent.size() == expected.max_independent && valid(g, independent, true));

    vector<int> maximal;
    g.maximal_cliques([&](const vector<int>& c){
        assert(valid(g, c, false));
        int mask = 0;
        for (int v : c) mask |= 1 << v;
        maximal.push_back(mask);
    });
    sort(maximal.begin(), maximal.end());
    assert(maximal == expected.maximal);
}

/// Larger graphs: every reported clique is maximal and distinct, and the largest equals max_clique
size_t largest_maximal(CliqueGraph& g){
    set<vector<int>> seen;
    size_t largest = 0;
    g.maximal_cliques([&](const vector<int>& c){
        assert(valid(g, c, false));
        CliqueGraph::Row common;
        for (int v = 0; v < g.n; v++) common[v] = 1;
        for (int v : c) common &= g.adj[v];
        assert(common.none());

        vector<int> key = c;
        sort(key.begin(), key.end());
        assert(seen.insert(key).second);
        largest = max(largest, c.size());
    });
    return largest;
}

CliqueGraph random_graph(int n, int density, vector<int>& adj){
    CliqueGraph g(n);
    adj.assign(n, 0);
    for (int u = 0; u < n; u++){
        for (int v = u + 1; v < n; v++){
            if (stress::rand_int(1, 100) > density) continue;
            g.add_edge(u, v);
            if (n <= 18) adj[u] |= 1 << v, adj[v] |= 1 << u;
        }
    }
    return g;
}

int main(){
    for (int n = 0; n <= 6; n++){
        vector<pair<int, int>> pairs;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) pairs.push_back({u, v});
        }

        for (int edges = 0; edges < (1 << pairs.size()); edges++){
            CliqueGraph g(n);
            vector<int> adj(n);
            for (int i = 0; i < (int)pairs.size(); i++){
                if (!(edges >> i & 1)) continue;
                auto [u, v] = pairs[i];
                g.add_edge(u, v), g.add_edge(v, u);
                adj[u] |= 1 << v, adj[v] |= 1 << u;
            }
            check(g, adj);
        }
    }

    for (long long it = 0; it < stress::scaled(1500); it++){
        int n = stress::rand_int(1, it % 8 ? 13 : 18), density = stress::rand_int(0, 100);
        vector<int> adj;
        CliqueGraph g = random_graph(n, density, adj);
        for (int v = 0; v < n; v++){
            if (stress::rand_int(0, 3) == 0) g.add_edge(v, v);
        }
        check(g, adj);

        for (int extra = stress::rand_int(0, 5); extra > 0; extra--){
            int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
            g.add_edge(u, v);
            if (u != v) adj[u] |= 1 << v, adj[v] |= 1 << u;
        }
        check(g, adj);
    }

    for (long long it = 0; it < stress::scaled(150); it++){
        /// Above density 0.75 n stops at 40 so Bron-Kerbosch, the reference here, stays fast
        int density = stress::rand_int(0, 95), n = stress::rand_int(19, density > 75 ? 40 : 50);
        vector<int> unused;
        CliqueGraph g = random_graph(n, density, unused);
        vector<int> clique = g.max_clique();
        assert(valid(g, clique, false) && clique.size() == largest_maximal(g));
    }

    /// n = 150 at density 0.7: the same answer size after relabelling the vertices
    {
        int n = 150;
        vector<int> unused, perm(n);
        CliqueGraph g = random_graph(n, 70, unused);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), stress::rng());

        CliqueGraph relabelled(n);
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (g.adj[u][v]) relabelled.add_edge(perm[u], perm[v]);
            }
        }
        vector<int> clique = g.max_clique();
        assert(valid(g, clique, false) && clique.size() == relabelled.max_clique().size());
    }

    /// n = CLIQUE_MAXN: a planted clique in a sparse graph, and a planted independent set in a dense one
    {
        int n = CLIQUE_MAXN;
        vector<int> unused, perm(n);
        iota(perm.begin(), perm.end(), 0);
        shuffle(perm.begin(), perm.end(), stress::rng());

        CliqueGraph sparse = random_graph(n, 10, unused);
        for (int i = 0; i < 20; i++){
            for (int j = 0; j < i; j++) sparse.add_edge(perm[i], perm[j]);
        }
        vector<int> clique = sparse.max_clique();
        assert(valid(sparse, clique, false) && clique.size() >= 20 && clique.size() == largest_maximal(sparse));

        CliqueGraph dense = random_graph(n, 90, unused);
        CliqueGraph complement(n);
        for (int i = 0; i < 20; i++){
            for (int j = 0; j < i; j++) dense.adj[perm[i]][perm[j]] = dense.adj[perm[j]][perm[i]] = 0;
        }
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++){
                if (!dense.adj[u][v]) complement.add_edge(u, v);
            }
        }
        vector<int> independent = dense.max_independent_set();
        assert(valid(dense, independent, true) && independent.size() >= 20 && independent.size() == largest_maximal(complement));
    }

    /// Labelled paths and cycles at n = CLIQUE_MAXN: their complements are ~99% dense with almost every degree equal,
    /// so a degree sort that scrambles ties leaves branch and bound exponential
    for (int cycle = 0; cycle < 2; cycle++){
        int n = CLIQUE_MAXN;
        CliqueGraph g(n);
        for (int v = 0; v + 1 < n; v++) g.add_edge(v, v + 1);
        if (cycle) g.add_edge(n - 1, 0);
        vector<int> independent = g.max_independent_set();
        assert(valid(g, independent, true) && (int)independent.size() == (cycle ? n / 2 : (n + 1) / 2));
    }
    return 0;
}
