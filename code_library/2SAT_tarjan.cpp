#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAXN 100010
#define MAXM 10000010

using namespace std;


struct Edge{
    int v, nxt;
} E[MAXM];

/***
 *
 * 2 SAT (1 based index for variables)
 * each variable can have two possible values - true or false
 * variables must satisfy a system of constraints on pairs of variables
 * positive variables are identified from 1 to n, while negative variables from -1 to -n
 *
***/

struct TwoSat{
    bool visited[MAXN * 2];
    int n, m, top, stk[MAXN * 2], head[MAXN * 2];

    TwoSat(int total_nodes){
        m = 0, top = 0;
        n = total_nodes * 2;
        memset(head, -1, sizeof(head));
    }

    inline void add_edge(int u, int v){
        E[m].v = v, E[m].nxt = head[u], head[u] = m++;
    }

    /// add implication, if a then b
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
        memset(visited, 0, sizeof(visited));

        for (int i = 2; i <= n; i += 2){
            if (!visited[i] && !visited[i ^ 1]) {
                top = 0;
                if (!dfs(i)){
                    while (top) visited[stk[--top]] = false;
                    if (!dfs(i ^ 1)) {
                        return false;
                    }
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
        for (int i = 2; i <= n * 2; i += 2){
            if (visited[i]) set_nodes.push_back(i >> 1);
        }

        return set_nodes;
    }
};


int main(){
    const int n = 5;
    auto sat = TwoSat(n);

    printf("%d\n\n", sat.is_satisfiable()); /// True

    sat.add_implication(-1, -2);            /// if 1 is false then 2 is false
    sat.add_implication(-2, -3);            /// if 2 is false then 3 is false
    sat.add_implication(-3, 4);             /// if 3 is false then 4 is true
    sat.add_implication(-3, -5);            /// if 3 is false then 5 is false
    sat.add_xor(1, 4);                      /// exactly one of 1 or 4 must be true
    sat.add_or(4, 5);                       /// either of 4 or 5 must be true

    printf("%d\n", sat.is_satisfiable());   /// True
    auto set_nodes = sat.get_assignment();
    for (auto x: set_nodes){
        printf("Node %d is set\n", x);      /// 1, 2, 3 and 5 will be true
    }
    puts("");

    sat.add_xor(1, 5);

    printf("%d\n", sat.is_satisfiable());   /// False
    return 0;
}
