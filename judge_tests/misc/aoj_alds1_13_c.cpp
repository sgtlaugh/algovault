// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_13_C
// competitive-verifier: TLE 3
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/15_puzzle_solver.cpp"
#undef main

int main(){
    vector<int> tiles(16);
    for (int& x : tiles){
        if (scanf("%d", &x) != 1) return 0;
    }

    FifteenPuzzle p(tiles);
    printf("%d\n", p.solve().first);
    return 0;
}
