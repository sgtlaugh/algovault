// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DPL_3_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/monotonic_stack.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<long long> h(n);
    for (long long& x: h){
        if (scanf("%lld", &x) != 1) return 0;
    }

    printf("%lld\n", largest_rectangle(h));
    return 0;
}
