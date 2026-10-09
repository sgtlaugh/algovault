#include "../common.h"

#define main library_main
#include "../../code_library/graphs/edge_coloring.cpp"
#undef main

/// KACTL EdgeColoring.h, unchanged apart from the macros it relies on
/// VizingEdgeColoring is a rewrite of it, so comparing against it is only an equivalence check: check() is the real oracle
namespace kactl{
    #define rep(i, a, b) for (int i = a; i < (b); ++i)
    #define all(x) begin(x), end(x)
    #define sz(x) (int)(x).size()
    typedef pair<int, int> pii;
    typedef vector<int> vi;

    vi edgeColoring(int N, vector<pii> eds) {
        vi cc(N + 1), ret(sz(eds)), fan(N), free(N), loc;
        for (pii e : eds) ++cc[e.first], ++cc[e.second];
        int u, v, ncols = *max_element(all(cc)) + 1;
        vector<vi> adj(N, vi(ncols, -1));
        for (pii e : eds) {
            tie(u, v) = e;
            fan[0] = v;
            loc.assign(ncols, 0);
            int at = u, end = u, d, c = free[u], ind = 0, i = 0;
            while (d = free[v], !loc[d] && (v = adj[u][d]) != -1)
                loc[d] = ++ind, cc[ind] = d, fan[ind] = v;
            cc[loc[d]] = c;
            for (int cd = d; at != -1; cd ^= c ^ d, at = adj[at][cd])
                swap(adj[at][cd], adj[end = at][cd ^ c ^ d]);
            while (adj[fan[i]][d] != -1) {
                int left = fan[i], right = fan[++i], e = cc[i];
                adj[u][e] = left;
                adj[left][e] = u;
                adj[right][e] = -1;
                free[right] = e;
            }
            adj[u][d] = fan[i];
            adj[fan[i]][d] = u;
            for (int y : {fan[0], u, end})
                for (int& z = free[y] = 0; adj[y][z] != -1; z++);
        }
        rep(i,0,sz(eds))
            for (tie(u, v) = eds[i]; adj[u][ret[i]] != v;) ++ret[i];
        return ret;
    }

    #undef rep
    #undef all
    #undef sz
}

/// Returns the number of distinct colors after asserting the coloring is proper and uses colors 0..limit - 1 only
int check(int n, const vector<array<int, 2>>& edges, const vector<int>& color, int limit){
    assert(color.size() == edges.size());
    vector<vector<int>> seen(n);
    set<int> distinct;
    for (int i = 0; i < (int)edges.size(); i++){
        assert(0 <= color[i] && color[i] < limit);
        distinct.insert(color[i]);
        for (int x : edges[i]) seen[x].push_back(color[i]);
    }
    for (auto& list : seen){
        sort(list.begin(), list.end());
        assert(adjacent_find(list.begin(), list.end()) == list.end());
    }
    return distinct.size();
}

int max_degree(int n, const vector<array<int, 2>>& edges){
    vector<int> deg(n);
    for (auto& e : edges) deg[e[0]]++, deg[e[1]]++;
    return n ? *max_element(deg.begin(), deg.end()) : 0;
}

/// Edges are given with right vertices already shifted by n_left
void run_bipartite(int n_left, int n_right, const vector<array<int, 2>>& edges){
    BipartiteEdgeColoring g(n_left, n_right);
    for (auto& e : edges) g.add_edge(e[0], e[1] - n_left);
    int d = max_degree(n_left + n_right, edges);
    assert(check(n_left + n_right, edges, g.solve(), d) == d);
}

vector<array<int, 2>> random_bipartite(int n_left, int n_right, int m){
    vector<array<int, 2>> edges;
    for (int i = 0; i < m && n_left && n_right; i++){
        edges.push_back({(int)stress::rand_int(0, n_left - 1), n_left + (int)stress::rand_int(0, n_right - 1)});
    }
    return edges;
}

void run_vizing(int n, const vector<array<int, 2>>& edges){
    VizingEdgeColoring g(n);
    vector<pair<int, int>> eds;
    for (auto& e : edges) g.add_edge(e[0], e[1]), eds.push_back({e[0], e[1]});
    vector<int> color = g.solve();
    check(n, edges, color, max_degree(n, edges) + 1);
    if (n) assert(color == kactl::edgeColoring(n, eds));
}

vector<array<int, 2>> random_simple(int n, int m){
    set<pair<int, int>> taken;
    vector<array<int, 2>> edges;
    m = min<long long>(m, 1LL * n * (n - 1) / 2);
    while ((int)edges.size() < m){
        int u = stress::rand_int(0, n - 1), v = stress::rand_int(0, n - 1);
        if (u == v || taken.count({min(u, v), max(u, v)})) continue;
        taken.insert({min(u, v), max(u, v)});
        edges.push_back({u, v});
    }
    return edges;
}

int main(){
    for (long long it = 0; it < stress::scaled(4000); it++){
        int n_left = stress::rand_int(0, 6), n_right = stress::rand_int(0, 6);
        run_bipartite(n_left, n_right, random_bipartite(n_left, n_right, stress::rand_int(0, 20)));
    }

    for (long long it = 0; it < stress::scaled(300); it++){
        int n_left = stress::rand_int(1, 4), n_right = stress::rand_int(1, 4);
        run_bipartite(n_left, n_right, random_bipartite(n_left, n_right, stress::rand_int(50, 300)));
    }

    for (long long it = 0; it < stress::scaled(60); it++){
        int n_left = stress::rand_int(1, 2000), n_right = stress::rand_int(1, 2000);
        run_bipartite(n_left, n_right, random_bipartite(n_left, n_right, stress::rand_int(0, 20000)));
    }

    /// K_{n,n}: every color class has to be a perfect matching
    for (int n = 1; n <= 120; n += 7){
        vector<array<int, 2>> edges;
        for (int u = 0; u < n; u++){
            for (int v = 0; v < n; v++) edges.push_back({u, n + v});
        }
        shuffle(edges.begin(), edges.end(), stress::rng());
        run_bipartite(n, n, edges);
    }

    for (long long it = 0; it < stress::scaled(4000); it++){
        int n = stress::rand_int(0, 8);
        run_vizing(n, random_simple(n, stress::rand_int(0, 28)));
    }

    for (long long it = 0; it < stress::scaled(40); it++){
        int n = stress::rand_int(2, 1500);
        run_vizing(n, random_simple(n, stress::rand_int(0, 30000)));
    }

    /// K_n for odd n has chromatic index n = D + 1, so every one of the D + 1 colors must appear
    for (int n = 1; n <= 45; n++){
        vector<array<int, 2>> edges;
        for (int u = 0; u < n; u++){
            for (int v = u + 1; v < n; v++) edges.push_back({u, v});
        }
        shuffle(edges.begin(), edges.end(), stress::rng());
        VizingEdgeColoring g(n);
        for (auto& e : edges) g.add_edge(e[0], e[1]);
        int used = check(n, edges, g.solve(), n);
        if (n % 2 == 1 && n > 1) assert(used == n);
    }

    return 0;
}
