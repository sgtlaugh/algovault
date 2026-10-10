// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/division_of_big_integers
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t;
    cin >> t;

    string out;
    while (t--){
        Bignum x, y;
        cin >> x >> y;
        out += (x / y).to_string();
        out += ' ';
        out += (x % y).to_string();
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
