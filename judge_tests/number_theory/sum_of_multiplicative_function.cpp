// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sum_of_multiplicative_function
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/min_25_sieve.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;
    while (t--){
        long long n, a, b;
        if (scanf("%lld %lld %lld", &n, &a, &b) != 3) return 0;

        Min25Sieve sieve(n, 469762049, {a, b});
        printf("%llu\n", (unsigned long long)sieve.sum([&](long long p, int e, long long){ return a * e + b * p; }));
    }
    return 0;
}
