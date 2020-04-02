#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAX 100010

using namespace std;


/***
 * 2 SAT (1 based index for variables)
 * Each variable can have two possible values - true or false
 * Variables must satisfy a system of constraints on pairs of variables
***/

namespace sat{
    bool visited[MAX * 2];
    vector <int> adj[MAX * 2], rev[MAX * 2];
    int n, m, l, dfs_t[MAX * 2], order[MAX * 2], parent[MAX * 2];

    inline int neg(int x){
        return ((x) <= n ? (x + n) : (x - n));
    }

    /// Call init once
    void init(int nodes){
        n = nodes, m = nodes * 2;
        for (int i = 0; i < MAX * 2; i++){
            adj[i].clear();
            rev[i].clear();
        }
    }

    /// Add implication, if a then b
    inline void add_implication(int a, int b){
        if (a < 0) a = n - a;
        if (b < 0) b = n - b;

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

    /// force variable x to be true (if x is negative, force !x to be true)
    inline void force_true(int x){
        if (x < 0) x = n - x;
        add_implication(neg(x), x);
    }

    /// force variable x to be false (if x is negative, force !x to be false)
    inline void force_false(int x){
        if (x < 0) x = n - x;
        add_implication(x, neg(x));
    }

    inline void topsort(int i){
        visited[i] = true;
        int j, x, len = rev[i].size();

        for (j = 0; j < len; j++){
            x = rev[i][j];
            if (!visited[x]) topsort(x);
        }
        dfs_t[i] = ++l;
    }

    inline void dfs(int i, int p){
        parent[i] = p;
        visited[i] = true;
        int j, x, len = adj[i].size();

        for (j = 0; j < len; j++){
            x = adj[i][j];
            if (!visited[x]) dfs(x, p);
        }
    }

    void build(){
        int i, x;
        memset(visited, 0, sizeof(visited));
        for (i = m, l = 0; i >= 1; i--){
            if (!visited[i]){
                topsort(i);
            }
            order[dfs_t[i]] = i;
        }

        memset(visited, 0, sizeof(visited));
        for (i = m; i >= 1; i--){
            x = order[i];
            if (!visited[x]){
                dfs(x, x);
            }
        }
    }

    /// Returns whether the system is 2-satisfiable
    bool is_satisfiable(){
        build();
        memset(visited, 0, sizeof(visited));

        for (int i = 1; i <= m; i++){
            int x = order[i];
            if (parent[x] == parent[neg(x)]) return false;
        }

        return true;
    }
}

int main(){
	int n = 4;
	sat::init(n);

	sat::add_implication(1, 2); 	/// if 1 is true then 2 is true
	sat::add_implication(-2, -3); 	/// if 2 is false then 3 is false
	sat::force_false(2); 			/// 2 must be false
	sat::add_xor(2, 4); 			/// exactly one of 2 or 4 must be true
	sat::add_or(1, 4); 				/// either 1 or 4 must be true

	puts(sat::is_satisfiable() ? "System is Satisfiable" : "System is not satisfiable");
    return 0;
}
