// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/two_sat
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/2SAT_lexicographic.cpp"
#undef main

int main(){
    int n, m;
    if (scanf(" p cnf %d %d", &n, &m) != 2) return 0;

    Graph g(n);
    for (int i = 0; i < m; i++){
        int a, b, z;
        if (scanf("%d %d %d", &a, &b, &z) != 3) return 0;
        g.add_or(a, b);
    }

    if (!g.is_satisfiable()){
        puts("s UNSATISFIABLE");
        return 0;
    }

    vector<char> truth(n + 1, 0);
    for (int x : g.get_assignment()) truth[x] = 1;

    string out = "s SATISFIABLE\nv";
    for (int x = 1; x <= n; x++) out += " " + to_string(truth[x] ? x : -x);
    out += " 0\n";
    fputs(out.c_str(), stdout);
    return 0;
}
