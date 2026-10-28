// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0085
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/josephus_problem.cpp"
#undef main

int main(){
    int n, m;
    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        printf("%lld\n", josephus(n, m, n));
    }

    return 0;
}
