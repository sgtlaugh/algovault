// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/primality_test
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/miller_rabin.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;
    while (q--){
        long long n;
        if (scanf("%lld", &n) != 1) return 0;
        puts(prm::is_prime(n) ? "Yes" : "No");
    }
    return 0;
}
