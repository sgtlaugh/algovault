// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/intersection_of_f2_vector_spaces
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/linear_algebra/xor_basis.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        XorBasis u, v;
        int n, m;
        ull x;
        if (scanf("%d", &n) != 1) return 0;
        for (int i = 0; i < n; i++){
            scanf("%llu", &x);
            u.insert(x);
        }
        if (scanf("%d", &m) != 1) return 0;
        for (int i = 0; i < m; i++){
            scanf("%llu", &x);
            v.insert(x);
        }

        XorBasis w = u.intersect(v);

        printf("%d", w.rank);
        for (ull b : w.basis){
            if (b) printf(" %llu", b);
        }
        putchar('\n');
    }
    return 0;
}
