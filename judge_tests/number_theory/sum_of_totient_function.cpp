// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sum_of_totient_function
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/du_sieve.cpp"
#undef main

int main(){
    long long n;
    if (scanf("%lld", &n) != 1) return 0;

    DuSieve du(n);
    printf("%lld\n", (long long)(du.phi_sum(n) % 998244353));
    return 0;
}
