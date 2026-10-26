// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/NTL_1_E
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/extended_gcd.cpp"
#undef main

int main(){
    long long a, b;
    if (scanf("%lld %lld", &a, &b) != 2) return 0;

    auto [g, x, y] = extended_gcd(a, b);
    printf("%lld %lld\n", x, y);
    return 0;
}
