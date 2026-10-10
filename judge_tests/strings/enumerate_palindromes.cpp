// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/enumerate_palindromes
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/manacher.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    auto pal = manacher(s);
    for (int i = 0; i < (int)pal.size(); i++) cout << pal[i] << (i + 1 < (int)pal.size() ? ' ' : '\n');
    return 0;
}
