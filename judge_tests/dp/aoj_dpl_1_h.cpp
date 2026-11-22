// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DPL_1_H
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/meet_in_the_middle.cpp"
#undef main

int main(){
    int n;
    long long capacity;
    if (scanf("%d %lld", &n, &capacity) != 2) return 0;

    vector<long long> values(n), weights(n);
    for (int i = 0; i < n; i++){
        if (scanf("%lld %lld", &values[i], &weights[i]) != 2) return 0;
    }

    printf("%lld\n", knapsack_large_weights(values, weights, capacity));
    return 0;
}
