/***
 *
 * Maximum Clique, Maximum Independent Set and Maximal Clique Enumeration
 * Finds a largest clique or independent set, or lists every maximal clique, of a graph with a few hundred vertices
 *
 * Complexity:
 *   - max_clique: exponential, at n = 150 on random graphs up to 2 s at edge density 0.9, 0.4 s at 0.8, 0.01 s at 0.5
 *   - max_independent_set: exponential, up to 0.15 s at n = 120 on paths and edge density 0.1, 16 s on paths at 160
 *   - maximal_cliques: O(3^(n/3) n^2 / 64), as a graph has at most 3^(n/3) maximal cliques
 *
 * CliqueGraph g(n): vertices 0..n - 1 for n <= 256, CliqueGraph<MAXN> g(n) takes n <= MAXN
 * g.add_edge(u, v): undirected edge, self-loops are ignored and duplicates are harmless
 * g.max_clique(): the vertices of one maximum clique, empty only for n = 0
 * g.max_independent_set(): the vertices of one maximum independent set, a maximum clique of the complement
 * g.maximal_cliques(callback): calls callback(const vector<int>& clique) once per maximal clique, in any vertex order
 *
 * Notes:
 *   - Every query can be called again, also after adding more edges
 *   - maximal_cliques reports an isolated vertex as a maximal clique, and the empty clique once for n = 0
 *   - Timings are measured with -O2
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template<int MAXN = 256>
struct CliqueGraph{
    using Row = bitset<MAXN>;

    struct Vertex{
        int id, bound;
    };

    /// Re-sort candidates by degree only while that work stays below this fraction of all expansions (MaxCliqueDyn)
    static constexpr double RESORT_LIMIT = 0.025;

    int n;
    vector<Row> adj;
    vector<vector<int>> color_class;
    vector<int> best, current;
    vector<long long> steps, old_steps;
    double total_steps;

    CliqueGraph(int n) : n(n), adj(n){
        assert(0 <= n && n <= MAXN);
    }

    void add_edge(int u, int v){
        if (u == v) return;
        adj[u][v] = adj[v][u] = 1;
    }

    vector<int> max_clique(){
        best.clear(), current.clear();
        if (n == 0) return best;

        color_class.assign(n + 1, {});
        steps.assign(n + 1, 0), old_steps.assign(n + 1, 0);
        total_steps = 0;

        vector<Vertex> all(n);
        for (int i = 0; i < n; i++) all[i] = {i, 0};
        sort_by_degree(all);
        expand(all, 1);
        return best;
    }

    vector<int> max_independent_set(){
        /// Bits at n and above are never read, so they need no masking
        CliqueGraph complement(n);
        for (int u = 0; u < n; u++){
            complement.adj[u] = ~adj[u];
            complement.adj[u][u] = 0;
        }
        return complement.max_clique();
    }

    template<typename F>
    void maximal_cliques(F callback){
        Row candidates, excluded;
        for (int i = 0; i < n; i++) candidates[i] = 1;

        vector<int> clique;
        bron_kerbosch(candidates, excluded, clique, callback);
    }

    void sort_by_degree(vector<Vertex>& r){
        for (auto& v : r){
            v.bound = 0;
            for (auto& u : r) v.bound += adj[v.id][u.id];
        }

        /// Stable: scrambled ties of near-regular graphs make the search exponential
        stable_sort(r.begin(), r.end(), [](const Vertex& a, const Vertex& b){ return a.bound > b.bound; });
        int max_degree = r[0].bound;
        for (int i = 0; i < (int)r.size(); i++) r[i].bound = min(i, max_degree) + 1;
    }

    void expand(vector<Vertex>& r, int level){
        steps[level] += steps[level - 1] - old_steps[level];
        old_steps[level] = steps[level - 1];

        while (!r.empty()){
            if (current.size() + r.back().bound <= best.size()) return;
            int v = r.back().id;
            current.push_back(v);

            vector<Vertex> next;
            for (auto& u : r){
                if (adj[v][u.id]) next.push_back({u.id, 0});
            }

            if (next.empty()){
                if (current.size() > best.size()) best = current;
            }
            else{
                if (steps[level]++ / ++total_steps < RESORT_LIMIT) sort_by_degree(next);
                color_sort(next);
                expand(next, level + 1);
            }
            current.pop_back(), r.pop_back();
        }
    }

    /// Greedy coloring: a vertex of color k can extend the clique by at most k, so colors below min_color are hopeless
    void color_sort(vector<Vertex>& r){
        int kept = 0, max_color = 1, min_color = max((int)best.size() - (int)current.size() + 1, 1);

        color_class[1].clear(), color_class[2].clear();
        for (auto& v : r){
            int k = 1;
            while (any_of(color_class[k].begin(), color_class[k].end(), [&](int u){ return adj[v.id][u]; })) k++;
            if (k > max_color) max_color = k, color_class[max_color + 1].clear();
            if (k < min_color) r[kept++].id = v.id;
            color_class[k].push_back(v.id);
        }

        /// The hopeless prefix ends in a zero bound, so expand stops once it pops down to it
        if (kept > 0) r[kept - 1].bound = 0;
        for (int k = min_color; k <= max_color; k++){
            for (int u : color_class[k]) r[kept++] = {u, k};
        }
    }

    template<typename F>
    void bron_kerbosch(Row candidates, Row excluded, vector<int>& clique, F& callback){
        if (candidates.none()){
            if (excluded.none()) callback(clique);
            return;
        }

        /// Tomita pivot: the vertex covering the most candidates, which is what yields the 3^(n/3) bound
        Row pool = candidates | excluded;
        int pivot = pool._Find_first();
        size_t covered = (candidates & adj[pivot]).count();
        for (int u = pool._Find_next(pivot); u < MAXN; u = pool._Find_next(u)){
            size_t count = (candidates & adj[u]).count();
            if (count > covered) pivot = u, covered = count;
        }

        Row branch = candidates & ~adj[pivot];
        for (int v = branch._Find_first(); v < MAXN; v = branch._Find_next(v)){
            clique.push_back(v);
            bron_kerbosch(candidates & adj[v], excluded & adj[v], clique, callback);
            clique.pop_back();
            candidates[v] = 0, excluded[v] = 1;
        }
    }
};

int main(){
    /// A house: square 0-1-2-3 with roof vertex 4 over edge 0-1, the roof the only triangle
    CliqueGraph house(5);
    for (auto [u, v] : vector<pair<int, int>>{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}}) house.add_edge(u, v);

    vector<int> roof = house.max_clique();
    sort(roof.begin(), roof.end());
    assert((roof == vector<int>{0, 1, 4}));
    assert(house.max_independent_set().size() == 2);  /// e.g. {2, 4}, no three vertices are pairwise apart

    vector<vector<int>> cliques;
    house.maximal_cliques([&](const vector<int>& c){
        cliques.push_back(c);
        sort(cliques.back().begin(), cliques.back().end());
    });
    sort(cliques.begin(), cliques.end());  /// every square edge but 0-1 is a maximal clique on its own
    assert((cliques == vector<vector<int>>{{0, 1, 4}, {0, 3}, {1, 2}, {2, 3}}));

    house.add_edge(0, 2);
    assert(house.max_clique().size() == 3);  /// now {0, 1, 2} ties the roof, no 4-clique
    return 0;
}
