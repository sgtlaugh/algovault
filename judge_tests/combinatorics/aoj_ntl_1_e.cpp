// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/NTL_1_E
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/combinatorics.cpp"
#undef main

int main(){
    long long a, b, x, y;
    if (scanf("%lld %lld", &a, &b) != 2) return 0;

    extended_gcd(a, b, x, y);
    printf("%lld %lld\n", x, y);
    return 0;
}
