// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/addition_of_big_integers
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
        out += (x + y).to_string();
        out += '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
