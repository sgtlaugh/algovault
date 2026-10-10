// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/stirling_number_of_the_second_kind
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/stirling_numbers.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    string out;
    char buf[16];
    for (long long x: stirling_second(n, 998244353)) out.append(buf, snprintf(buf, sizeof(buf), "%lld ", x));
    out.back() = '\n';

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
