// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/combinatorics.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != EOF && !isdigit(c)) c = getchar();

    long long x = 0;
    for (; isdigit(c); c = getchar()) x = x * 10 + (c - '0');
    return x;
}

int main(){
    int t = read_int();
    long long m = read_int();
    Combinatorics comb(min(m, 10000000LL) - 1, m);

    string out;
    char buf[16];
    while (t--){
        int n = read_int(), k = read_int();
        out.append(buf, snprintf(buf, sizeof(buf), "%lld\n", comb.nCr(n, k)));
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
