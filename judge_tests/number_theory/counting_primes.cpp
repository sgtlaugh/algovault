// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/counting_primes
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/fast_prime_counting.cpp"
#undef main

int main(){
    long long n;
    if (scanf("%lld", &n) != 1) return 0;

    PrimeCounter<> pc;
    printf("%llu\n", (unsigned long long)pc.lehmer(n));
    return 0;
}
