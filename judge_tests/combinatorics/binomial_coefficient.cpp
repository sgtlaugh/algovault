// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/binomial_coefficient
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/binomial_coefficients.cpp"
#undef main

int main(){
    int t;
    long long m;
    if (scanf("%d %lld", &t, &m) != 2) return 0;

    Binomial binomial(m);
    while (t--){
        long long n, k;
        if (scanf("%lld %lld", &n, &k) != 2) return 0;
        printf("%lld\n", (long long)binomial.binomial(n, k));
    }
    return 0;
}
