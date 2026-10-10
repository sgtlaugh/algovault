// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0085
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/josephus_problem.cpp"
#undef main

int main(){
    int n, m;
    while (scanf("%d %d", &n, &m) == 2 && (n || m)){
        int survivor = josephus1(n, m, n);
        assert(survivor == josephus2(n, m, n));
        printf("%d\n", survivor);
    }
    return 0;
}
