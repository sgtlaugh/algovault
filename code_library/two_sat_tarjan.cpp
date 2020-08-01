/***
 *
 * 2SAT with Tarjan's algorithm (1 based index for variables)
 * Each variable can have two possible values: true or false
 * Variables must satisfy a system of constraints on pairs of variables
 * Positive variables are identified from 1 to n, while negative variables from -1 to -n
 *
 * Complexity: O(n + m), n = number of nodes and m = number of edges or constraints
 *
***/



#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAXM 10000010

using namespace std;

struct Edge{
    int v, nxt;
} E[MAXM];


struct Graph{
    int n, m, top, lim;
    vector <int> stk, head, visited;

    Graph(int n) : n(2 * n) {
        m = 0, lim = this->n + 2;
        head.resize(lim, -1), visited.resize(lim, 0), stk.resize(lim, 0);
    }

    inline void add_edge(int u, int v){
        E[m].v = v, E[m].nxt = head[u], head[u] = m++;
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        assert(abs(a) != abs(b));

        a *= 2, b *= 2;
        if (a < 0) a = -a ^ 1;
        if (b < 0) b = -b ^ 1;
        add_edge(a, b);
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

    inline bool dfs(int u){
        if (visited[u ^ 1]) return false;
        if (visited[u]) return true;

        visited[u] = true, stk[top++] = u;
        for (int i = head[u]; i >= 0; i = E[i].nxt) {
            if (!dfs(E[i].v)) return false;
        }

        return true;
    }

    /// returns true if the system is 2-satisfiable and false otherwise
    bool is_satisfiable(){
        for (int i = 0; i < lim; i++) visited[i] = 0;

        for (int i = 2; i <= n; i += 2){
            if (!visited[i] && !visited[i ^ 1]) {
                top = 0;
                if (!dfs(i)){
                    while (top) visited[stk[--top]] = 0;
                    if (!dfs(i ^ 1)) return false;
                }
            }
        }

        return true;
    }

    /***
     * returns nodes set to true in a vector
     * only call when is_satisfiable()
    ***/

    vector <int> get_assignment(){
        vector <int> set_nodes;
        for (int i = 2; i <= n; i += 2){
            if (visited[i]) set_nodes.push_back(i >> 1);
        }

        return set_nodes;
    }
};


int main(){
    auto g = Graph(5);
    assert(g.is_satisfiable());

    g.add_implication(-1, -2);  /// if 1 is false then 2 is false
    g.add_implication(-2, -3);  /// if 2 is false then 3 is false
    g.add_implication(-3, 4);   /// if 3 is false then 4 is true
    g.add_implication(-3, -5);  /// if 3 is false then 5 is false
    g.add_xor(1, 4);            /// exactly one of 1 or 4 must be true
    g.add_or(4, 5);             /// either of 4 or 5 must be true

    assert(g.is_satisfiable());
    auto set_nodes = g.get_assignment();
    vector <int> expected_set_nodes{1, 2, 3, 5};
    assert(set_nodes == expected_set_nodes);

    g.add_xor(1, 5);
    assert(!g.is_satisfiable());
    return 0;
}
