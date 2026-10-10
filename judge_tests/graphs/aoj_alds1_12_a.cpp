// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_12_A
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/minimum_spanning_tree.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<vector<int>> w(n, vector<int>(n));
    for (auto& row : w){
        for (auto& x : row){
            if (scanf("%d", &x) != 1) return 0;
            if (x == -1) x = numeric_limits<int>::max();
        }
    }

    printf("%lld\n", dense_prim(w).weight);
    return 0;
}
