// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_tree_vertex_add_path_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/link_cut_tree.cpp"
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
    LinkCutTree<long long> lct(n);
    vector<long long> val(n);
    for (int i = 0; i < n; i++){
        val[i] = read_int();
        lct.set(i, val[i]);
    }

    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        lct.link(u, v);
    }

    while (q--){
        int t = read_int();
        if (t == 0){
            int u = read_int(), v = read_int(), w = read_int(), x = read_int();
            lct.cut(u, v);
            lct.link(w, x);
        }
        else if (t == 1){
            int p = read_int(), x = read_int();
            lct.set(p, val[p] += x);
        }
        else{
            int u = read_int(), v = read_int();
            printf("%lld\n", lct.query(u, v));
        }
    }
    return 0;
}
