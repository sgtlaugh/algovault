// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_10_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/hirschberg_lcs.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    static char x[1010], y[1010];
    while (q--){
        if (scanf("%1009s %1009s", x, y) != 2) break;
        printf("%zu\n", lcs_sequence(string(x), string(y)).size());
    }

    return 0;
}
