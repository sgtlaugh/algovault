// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/wildcard_pattern_matching
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/wildcard_matching.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s, t;
    cin >> s >> t;

    string res(s.size() - t.size() + 1, '0');
    for (int i : wildcard_match(s, t, '*')) res[i] = '1';
    cout << res << '\n';
    return 0;
}
