// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2286
// competitive-verifier: TLE 8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/mobius_function.cpp"
#undef main

/// F_n holds 0/1 and one fraction a/b per coprime pair 1 <= a <= b <= n; ordered coprime pairs number sum mu(d) (n/d)^2
int main(){
    const int N = 1000000;
    MobiusSieve sieve(N);
    vector<int> prefix(N + 1, 0);
    for (int d = 1; d <= N; d++) prefix[d] = prefix[d - 1] + sieve.mu[d];

    int t, n;
    if (scanf("%d", &t) != 1) return 0;
    while (t-- && scanf("%d", &n) == 1){
        long long pairs = 0;
        for (int l = 1, r; l <= n; l = r + 1){
            long long q = n / l;
            r = n / q;
            pairs += (prefix[r] - prefix[l - 1]) * q * q;
        }
        printf("%lld\n", 1 + (pairs + 1) / 2);
    }

    return 0;
}
