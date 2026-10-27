// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/factorial
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/fast_factorial.cpp"
#undef main

int main(){
    int t;
    cin >> t;

    FastFactorial fact(998244353, t);
    while (t--){
        long long n;
        cin >> n;
        cout << fact(n) << '\n';
    }

    return 0;
}
