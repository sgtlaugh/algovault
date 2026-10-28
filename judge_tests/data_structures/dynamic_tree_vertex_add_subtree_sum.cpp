// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_tree_vertex_add_subtree_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/euler_tour_tree.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<long long> a(n);
    for (auto& x : a) x = read_int();

    EulerTourTree<long long> ett(a);
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        ett.link(u, v);
    }

    while (q--){
        int t = read_int();
        if (t == 0){
            int u = read_int(), v = read_int(), w = read_int(), x = read_int();
            ett.cut(u, v);
            ett.link(w, x);
        }
        else if (t == 1){
            int p = read_int(), x = read_int();
            ett.add(p, x);
        }
        else{
            int v = read_int(), p = read_int();
            printf("%lld\n", ett.subtree_query(v, p));
        }
    }

    return 0;
}
