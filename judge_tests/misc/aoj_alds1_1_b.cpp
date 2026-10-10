// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_1_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/assembly.cpp"
#undef main

int main(){
    int a, b;
    if (scanf("%d %d", &a, &b) != 2) return 0;

    printf("%u\n", asm_gcd(a, b));
    return 0;
}
