// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/suffixarray
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/suffix_array.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    auto sa = suffix_array(s).sa;
    for (int i = 0; i < (int)sa.size(); i++) cout << sa[i] << (i + 1 < (int)sa.size() ? ' ' : '\n');
    return 0;
}
