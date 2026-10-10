// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lyndon_factorization
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/lyndon_factorization.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    for (int i : lyndon_factorization(s)) cout << i << ' ';
    cout << s.size() << '\n';
    return 0;
}
