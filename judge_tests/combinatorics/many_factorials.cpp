// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/many_factorials
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/fast_factorial.cpp"
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
    FastFactorial fact(998244353, t);

    string out;
    char buf[16];
    while (t--){
        long long n = read_int();
        out.append(buf, snprintf(buf, sizeof(buf), "%lld\n", fact(n)));
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
