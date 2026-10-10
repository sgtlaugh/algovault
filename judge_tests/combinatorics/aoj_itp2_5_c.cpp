// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ITP2_5_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/permutation_rank.cpp"
#undef main

void print(const vector<int>& a){
    for (int i = 0; i < (int)a.size(); i++) printf("%d%c", a[i], i + 1 == (int)a.size() ? '\n' : ' ');
}

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<int> a(n);
    for (auto& x: a){
        if (scanf("%d", &x) != 1) return 0;
    }

    long long rank = find_rank(a), total = 1;
    for (int i = 2; i <= n; i++) total *= i;

    if (rank > 1) print(find_permutation(n, rank - 1));
    print(a);
    if (rank < total) print(find_permutation(n, rank + 1));
    return 0;
}
