// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/dynamic_tree_vertex_set_path_composite
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/link_cut_tree.cpp"
#undef main

const long long MOD = 998244353;

/// f(x) = a x + b along the path in order, with the reverse order kept for re-rooting
struct Fn{
    long long a, b, ra, rb;
};

struct Compose{
    static Fn combine(Fn f, Fn g){
        return {g.a * f.a % MOD, (g.a * f.b + g.b) % MOD, f.ra * g.ra % MOD, (f.ra * g.rb + f.rb) % MOD};
    }
    static Fn apply(Fn val, int, int){ return val; }
    static int compose(int, int lz){ return lz; }
    static Fn reversed(Fn f){ return {f.ra, f.rb, f.a, f.b}; }
};

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    LinkCutTree<Fn, Compose, int> lct(n);
    for (int i = 0; i < n; i++){
        long long a = read_int(), b = read_int();
        lct.set(i, {a, b, a, b});
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
            long long p = read_int(), c = read_int(), d = read_int();
            lct.set(p, {c, d, c, d});
        }
        else{
            int u = read_int(), v = read_int();
            long long x = read_int();
            Fn f = lct.query(u, v);
            printf("%lld\n", (f.a * x + f.b) % MOD);
        }
    }

    return 0;
}
