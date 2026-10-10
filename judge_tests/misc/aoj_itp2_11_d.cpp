// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ITP2_11_D
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/bit_twiddling.cpp"
#undef main

int main(){
    int n, k;
    if (scanf("%d %d", &n, &k) != 2) return 0;

    /// next_same_popcount needs a nonzero word, the only 0-subset is printed directly
    if (k == 0){
        puts("0:");
        return 0;
    }

    for (unsigned int x = (1u << k) - 1; x < (1u << n); x = next_same_popcount(x)){
        printf("%u:", x);
        for (int i : set_bit_indices(x)) printf(" %d", i);
        putchar('\n');
    }
    return 0;
}
