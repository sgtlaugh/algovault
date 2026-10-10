// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/NTL_2_E
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

int main(){
    Bignum x, y;
    cin >> x >> y;
    cout << x % y << '\n';
    return 0;
}
