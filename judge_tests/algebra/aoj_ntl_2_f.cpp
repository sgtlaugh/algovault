// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/NTL_2_F
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    Bignum a, b;
    cin >> a >> b;
    cout << a * b << '\n';
    return 0;
}
