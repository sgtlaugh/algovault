// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/kth_root_integer
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/integer_root.cpp"
#undef main

unsigned long long read_uint(){
    int c = getchar();
    while (c < '0') c = getchar();
    unsigned long long x = 0;
    for (; c >= '0'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int t = read_uint();
    while (t--){
        unsigned long long a = read_uint();
        int k = read_uint();

        unsigned long long res = k == 2 ? isqrt(a) : k == 3 ? icbrt(a) : iroot(a, k);
        printf("%llu\n", res);
    }
    return 0;
}
