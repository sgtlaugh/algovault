/***
 *
 * 2SAT with Kosaraju's algorithm (1 based index for variables)
 * Each variable can have two possible values: true or false
 * Variables must satisfy a system of constraints on pairs of variables
 *
 * Complexity: O(n + m), n = number of nodes and m = number of edges or constraints
 *
***/


#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

struct Graph{
    int n, l;
    vector <vector<int>> adj, rev;
    vector <int> visited, parent, order, dfs_t;

    Graph() {}
    Graph(int n): n(n){
        int m = 2 * n + 2;
        adj.resize(m, vector<int>()), rev.resize(m, vector<int>());
        visited.resize(m, 0), dfs_t.resize(m, 0), order.resize(m, 0), parent.resize(m, 0);
    }

    inline int neg(int x){
        return ((x) <= n ? (x + n) : (x - n));
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        if (a < 0) a = n - a;
        if (b < 0) b = n - b;
        assert(a >= 1 && a <= 2 * n && b >= 1 && b <= 2 * n);
        adj[a].push_back(b);
        rev[b].push_back(a);
    }

    inline void add_or(int a, int b){
        add_implication(-a, b);
        add_implication(-b, a);
    }

    inline void add_xor(int a, int b){
        add_or(a, b);
        add_or(-a, -b);
    }

    inline void add_and(int a, int b){
        add_or(a, b);
        add_or(a, -b);
        add_or(-a, b);
    }

    /// Force variable x to be true (if x is negative, force !x to be true)
    inline void force_true(int x){
        if (x < 0) x = n - x;
        add_implication(neg(x), x);
    }

    /// Force variable x to be false (if x is negative, force !x to be false)
    inline void force_false(int x){
        if (x < 0) x = n - x;
        add_implication(x, neg(x));
    }

    inline void topsort(int i){
        visited[i] = true;
        for (auto x: rev[i]){
            if (!visited[x]) topsort(x);
        }
        dfs_t[i] = ++l;
    }

    inline void dfs(int i, int p){
        parent[i] = p, visited[i] = true;
        for (auto x: adj[i]){
            if (!visited[x]) dfs(x, p);
        }
    }

    void build(){
        int i, x;
        for (i = 0; i <= 2 * n; i++) visited[i] = 0;
        for (i = 2 * n, l = 0; i >= 1; i--){
            if (!visited[i]) topsort(i);
            order[dfs_t[i]] = i;
        }

        for (i = 0; i <= 2 * n; i++) visited[i] = 0;
        for (i = 2 * n; i >= 1; i--){
            x = order[i];
            if (!visited[x]) dfs(x, x);
        }
    }

    /// Returns whether the system is 2-satisfiable
    bool is_satisfiable(){
        build();
        for (int i = 0; i <= 2 * n; i++) visited[i] = 0;

        for (int i = 1; i <= 2 * n; i++){
            int x = order[i];
            if (parent[x] == parent[neg(x)]) return false;
        }

        return true;
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
    return 0;
}
