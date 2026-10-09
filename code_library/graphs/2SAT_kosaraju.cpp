/***
 *
 * 2SAT with Kosaraju's algorithm (1 based index for variables)
 * Each variable can have two possible values: true or false
 * Variables must satisfy a system of constraints on pairs of variables
 * After is_satisfiable() returns true, value(x) gives a satisfying assignment
 *
 * Complexity: O(n + m), n = number of nodes and m = number of edges or constraints
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct Graph{
    int n, t;
    vector<pair<int, int>> edges;
    vector<int> adj_start, adj_to, rev_start, rev_to;
    vector<int> visited, comp, order, dfs_t;

    Graph(int n): n(n){
        int m = 2 * n + 2;
        visited.resize(m, 0), dfs_t.resize(m, 0), order.resize(m, 0), comp.resize(m, 0);
    }

    inline int neg(int x){
        return ((x) <= n ? (x + n) : (x - n));
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        if (a < 0) a = n - a;
        if (b < 0) b = n - b;
        assert(a >= 1 && a <= 2 * n && b >= 1 && b <= 2 * n);
        edges.push_back({a, b}), edges.push_back({neg(b), neg(a)});
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    inline void add_and(int a, int b){
        force_true(a);
        force_true(b);
    }

    /// Force variable x to be true (if x is negative, force !x to be true)
    inline void force_true(int x){
        if (x < 0) x = n - x;
        edges.push_back({neg(x), x});  /// !x -> x is its own contrapositive, add it once
    }

    /// Force variable x to be false (if x is negative, force !x to be false)
    inline void force_false(int x){
        if (x < 0) x = n - x;
        edges.push_back({x, neg(x)});  /// x -> !x is its own contrapositive, add it once
    }

    /// Stable counting sort of the edges by source, so the DFS visits neighbours in insertion order
    void flatten(vector<int>& start, vector<int>& to, bool reverse){
        start.assign(2 * n + 2, 0), to.resize(edges.size());
        for (auto [u, v]: edges) start[reverse ? v : u]++;
        for (int i = 1; i <= 2 * n + 1; i++) start[i] += start[i - 1];
        for (int i = (int)edges.size() - 1; i >= 0; i--){
            auto [u, v] = edges[i];
            if (reverse) swap(u, v);
            to[--start[u]] = v;
        }
    }

    inline void topsort(int i){
        visited[i] = true;
        for (int e = rev_start[i]; e < rev_start[i + 1]; e++){
            if (!visited[rev_to[e]]) topsort(rev_to[e]);
        }
        dfs_t[i] = ++t;
    }

    inline void dfs(int i, int c){
        comp[i] = c, visited[i] = true;
        for (int e = adj_start[i]; e < adj_start[i + 1]; e++){
            if (!visited[adj_to[e]]) dfs(adj_to[e], c);
        }
    }

    /// Components are numbered in discovery order, which is a reverse topological order of the implication graph
    void build(){
        int i, x, c = 0;
        flatten(adj_start, adj_to, false), flatten(rev_start, rev_to, true);

        for (i = 0; i <= 2 * n; i++) visited[i] = 0;
        for (i = 2 * n, t = 0; i >= 1; i--){
            if (!visited[i]) topsort(i);
            order[dfs_t[i]] = i;
        }

        for (i = 0; i <= 2 * n; i++) visited[i] = 0;
        for (i = 2 * n; i >= 1; i--){
            x = order[i];
            if (!visited[x]) dfs(x, c++);
        }
    }

    /// Returns whether the system is 2-satisfiable
    bool is_satisfiable(){
        build();
        for (int i = 1; i <= n; i++){
            if (comp[i] == comp[i + n]) return false;
        }
        return true;
    }

    /// Value of literal x in a satisfying assignment, valid only after is_satisfiable() returned true
    /// If x implies !x then comp[!x] is discovered first, so x is false
    inline bool value(int x){
        if (x < 0) return !value(-x);
        return comp[x] < comp[neg(x)];
    }
};


int main(){
    auto g = Graph(4);

    g.add_implication(1, 2);     /// if 1 is true then 2 is true
    g.add_implication(-2, -3);   /// if 2 is false then 3 is false
    g.force_false(2);            /// 2 must be false
    g.add_xor(2, 4);             /// exactly one of 2 or 4 must be true
    g.add_or(1, 4);              /// either 1 or 4 must be true

    assert(g.is_satisfiable());
    assert(!g.value(1) && !g.value(2) && !g.value(3) && g.value(4) && !g.value(-4));
    return 0;
}
