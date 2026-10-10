// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_10_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/hunt_szymanski.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int q;
    cin >> q;
    string x, y;
    while (q--){
        cin >> x >> y;
        cout << lcs(x.c_str(), y.c_str()) << '\n';
    }
    return 0;
}
