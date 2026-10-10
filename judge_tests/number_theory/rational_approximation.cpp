// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/rational_approximation
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/continued_fractions.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long n, x, y;
        if (scanf("%lld %lld %lld", &n, &x, &y) != 3) return 0;

        auto upper = fraction_search([&](long long p, long long q){ return p * y >= x * q; }, n);
        /// The largest p / q <= x / y is the smallest q / p >= y / x, so search the reciprocal and flip it back
        auto lower = fraction_search([&](long long q, long long p){ return q * x >= y * p; }, n);
        printf("%lld %lld %lld %lld\n", lower.second, lower.first, upper.first, upper.second);
    }
    return 0;
}
