// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/zalgorithm
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/z_algorithm.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    auto z = z_function(s);
    for (int i = 0; i < (int)z.size(); i++) cout << z[i] << (i + 1 < (int)z.size() ? ' ' : '\n');
    return 0;
}
