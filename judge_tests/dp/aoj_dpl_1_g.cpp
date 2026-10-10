// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DPL_1_G
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/bounded_knapsack.cpp"
#undef main

int main(){
    int n, capacity;
    if (scanf("%d %d", &n, &capacity) != 2) return 0;

    vector<long long> values(n), counts(n);
    vector<int> weights(n);
    for (int i = 0; i < n; i++){
        if (scanf("%lld %d %lld", &values[i], &weights[i], &counts[i]) != 3) return 0;
    }

    printf("%lld\n", bounded_knapsack(values, weights, counts, capacity));
    return 0;
}
