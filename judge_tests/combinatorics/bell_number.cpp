// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/bell_number
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/bell_numbers.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    string out;
    char buf[16];
    for (long long x: bell_numbers_ntt(n)) out.append(buf, snprintf(buf, sizeof(buf), "%lld ", x));
    out.back() = '\n';

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
