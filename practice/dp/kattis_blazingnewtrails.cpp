/***
 *
 * Aliens Trick (WQS Binary Search, Lagrangian Relaxation)
 * Minimum cost of a solution with exactly k items, when that cost f(k) is convex in k
 *
 * Complexity: O(log(hi - lo)) calls to the solver, plus one
 *
 * aliens_trick(k, lo, hi, solver) returns f(k)
 *   solver(lambda) returns {value, count}: value = min over all solutions of cost + lambda * items,
 *   count = the MINIMUM number of items among the solutions reaching that value
 *   (a DP over pairs {cost, items} compared lexicographically gives exactly this)
 *
 * Requirements:
 *   - integer costs, f convex on its domain and k inside the domain
 *   - [lo, hi] contains every marginal -(f(c) - f(c - 1)), so M = max |f(c) - f(c - 1)| gives [-M, M]
 *   - the bounds' type is wide enough for hi - lo, the solver's value type for cost + max(|lo|, |hi|) * (max items)
 *   - f(k) comes back in the solver's value type, so int bounds with an int64 solver still return int64
 *
 * Ties: several counts can be optimal for one lambda, so a lambda with count == k may not exist
 * The driver takes the smallest lambda whose minimum optimal count is <= k; k is then optimal there too,
 * and f(k) = value - lambda * k (never value - lambda * count)
 *
 * Maximizing a concave f: minimize the negated gains and negate the result
 *
 * Example: split a non-negative array into k segments minimizing the sum of squared segment sums
 *   aliens_trick(k, 0LL, S * S, solver), solver(lambda) runs the O(n^2) split DP with + lambda per segment
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Solver>
auto aliens_trick(int k, T lo, T hi, Solver solver){
    while (lo < hi){
        T mid = lo + (hi - lo) / 2;
        if (solver(mid).second <= k) hi = mid;
        else lo = mid + 1;
    }

    auto value = solver(lo).first;
    return value - decltype(value)(lo) * k;
}

struct DSU{
    vector<int> par;

    explicit DSU(int n) : par(n){
        iota(par.begin(), par.end(), 0);
    }

    int find(int v){
        while (v != par[v]) v = par[v] = par[par[v]];
        return v;
    }

    bool unite(int a, int b){
        a = find(a), b = find(b);
        if (a == b) return false;
        par[a] = b;
        return true;
    }
};

struct Edge{
    int u, v, c;
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m, k, w;
    cin >> n >> m >> k >> w;
    vector<bool> special(n + 1, false);
    for (int i = 0; i < k; i++){
        int s;
        cin >> s;
        special[s] = true;
    }

    vector<Edge> cross, plain;
    for (int i = 0; i < m; i++){
        Edge e;
        cin >> e.u >> e.v >> e.c;
        (special[e.u] != special[e.v] ? cross : plain).push_back(e);
    }
    auto by_cost = [](const Edge& a, const Edge& b){ return a.c < b.c; };
    sort(cross.begin(), cross.end(), by_cost);
    sort(plain.begin(), plain.end(), by_cost);

    /// Kruskal with lambda added to every cross edge; both lists stay sorted under the shift, so a merge
    /// replaces the sort, and taking plain edges first on ties yields the fewest cross edges among the MSTs
    int merged = 0;
    auto solver = [&](long long lambda){
        DSU dsu(n + 1);
        pair<long long, int> res = {0, 0};
        merged = 0;
        size_t i = 0, j = 0;
        while (i < cross.size() || j < plain.size()){
            bool take_cross = j == plain.size() || (i < cross.size() && cross[i].c + lambda < plain[j].c);
            const Edge& e = take_cross ? cross[i++] : plain[j++];
            if (!dsu.unite(e.u, e.v)) continue;
            merged++;
            res.first += e.c + (take_cross ? lambda : 0);
            res.second += take_cross;
        }
        return res;
    };

    /// Swapping one tree edge moves the cost by under 10^5, so these lambdas force the fewest and the most cross edges
    const long long BOUND = 100000;
    int fewest = solver(BOUND).second;
    if (merged != n - 1){
        cout << -1 << '\n';
        return 0;
    }
    int most = solver(-BOUND).second;
    if (w < fewest || w > most) cout << -1 << '\n';
    else cout << aliens_trick(w, -BOUND, BOUND, solver) << '\n';
    return 0;
}
