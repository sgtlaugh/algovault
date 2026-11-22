// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2286
// competitive-verifier: TLE 8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/linear_sieve.cpp"
#undef main

int main(){
    const int N = 1000000;
    LinearSieve sieve(N, LinearSieve::PHI);
    vector<long long> terms(N + 1, 1);
    for (int k = 1; k <= N; k++) terms[k] = terms[k - 1] + sieve.phi(k);

    int t, n;
    if (scanf("%d", &t) != 1) return 0;
    while (t-- && scanf("%d", &n) == 1) printf("%lld\n", terms[n]);

    return 0;
}
