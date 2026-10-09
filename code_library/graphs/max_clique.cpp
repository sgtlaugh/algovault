/***
 *
 * Maximum Clique, Maximum Independent Set and Maximal Clique Enumeration
 * Branch and bound over bitset adjacency (MaxCliqueDyn) and Bron-Kerbosch with Tomita pivoting
 *
 * Complexity: exponential worst case for all three
 *     max_clique: measured at n = 150 with -O2 on random graphs: 0.7-2 s at edge density 0.9 (the hard case),
 *         0.2-0.4 s at 0.8, under 0.01 s at 0.5
 *     max_independent_set: a sparse graph has a 90-99% dense complement, so keep sparse inputs near n <= 120
 *         measured with -O2: n = 120 shuffled paths and G(n, p <= 0.1) under 0.15 s, n = 160 shuffled paths 4-16 s
 *     maximal_cliques: O(3^(n/3) * n^2 / 64), a graph has at most 3^(n/3) maximal cliques (Moon-Moser)
 *
 * CliqueGraph g(n): vertices 0..n - 1, n <= CLIQUE_MAXN (raise the constant for bigger graphs)
 * g.add_edge(u, v): undirected edge, self-loops are ignored and duplicates are harmless
 * g.max_clique(): the vertices of one maximum clique, empty only for n = 0
 * g.max_independent_set(): the vertices of one maximum independent set, a maximum clique of the complement
 * g.maximal_cliques(callback): calls callback(const vector<int>& clique) once per maximal clique,
 *     vertices in no particular order, an isolated vertex is a maximal clique, n = 0 reports the empty clique once
 *
 * Every query can be called again, also after adding more edges
 *
 * Example:
 *     CliqueGraph g(4);
 *     g.add_edge(0, 1), g.add_edge(1, 2), g.add_edge(0, 2), g.add_edge(2, 3);
 *     g.max_clique();  // {0, 1, 2} in some order
 *     g.maximal_cliques([](const vector<int>& c){ ... });  // {0, 1, 2} and {2, 3}
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

const int CLIQUE_MAXN = 256;

struct CliqueGraph{
    using Row = bitset<CLIQUE_MAXN>;

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
        assert(0 <= n && n <= CLIQUE_MAXN);
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
        for (int u = pool._Find_next(pivot); u < CLIQUE_MAXN; u = pool._Find_next(u)){
            size_t count = (candidates & adj[u]).count();
            if (count > covered) pivot = u, covered = count;
        }

        Row branch = candidates & ~adj[pivot];
        for (int v = branch._Find_first(); v < CLIQUE_MAXN; v = branch._Find_next(v)){
            clique.push_back(v);
            bron_kerbosch(candidates & adj[v], excluded & adj[v], clique, callback);
            clique.pop_back();
            candidates[v] = 0, excluded[v] = 1;
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

    void sort_by_degree(vector<Vertex>& r){
        for (auto& v : r){
            v.bound = 0;
            for (auto& u : r) v.bound += adj[v.id][u.id];
        }

        /// Stable so near-regular graphs keep their labelling: scrambled ties made 256-vertex labelled paths exponential
        stable_sort(r.begin(), r.end(), [](const Vertex& a, const Vertex& b){ return a.bound > b.bound; });
        int max_degree = r[0].bound;
        for (int i = 0; i < (int)r.size(); i++) r[i].bound = min(i, max_degree) + 1;
    }
};

int main(){
    /***
     * House: square 0-1-2-3 with roof vertex 4 on top of edge 0-1
     * The roof triangle is the only triangle, every square edge but 0-1 is a maximal clique on its own
    ***/
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
    sort(cliques.begin(), cliques.end());
    assert((cliques == vector<vector<int>>{{0, 1, 4}, {0, 3}, {1, 2}, {2, 3}}));

    house.add_edge(0, 2);
    assert(house.max_clique().size() == 3);  /// now {0, 1, 2} ties the roof, no 4-clique
    return 0;
}
