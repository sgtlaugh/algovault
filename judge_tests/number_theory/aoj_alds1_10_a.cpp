// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_10_A
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/fast_fibonacci.cpp"
#undef main

/// AOJ starts the sequence at fib(0) = fib(1) = 1, one index ahead of F(0) = 0, F(1) = 1
int main(){
    long long n;
    if (scanf("%lld", &n) != 1) return 0;

    long long res = fibonacci(n + 1);
    assert(fibonacci(n + 1, LLONG_MAX) == res);
    printf("%lld\n", res);
    return 0;
}
