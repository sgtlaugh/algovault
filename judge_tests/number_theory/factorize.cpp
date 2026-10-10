// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/factorize
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/pollard_rho.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    PollardRho f;
    while (q--){
        long long a;
        if (scanf("%lld", &a) != 1) return 0;
        auto factors = f.factorize(a);
        printf("%d", (int)factors.size());
        for (long long p : factors) printf(" %lld", p);
        puts("");
    }
    return 0;
}
