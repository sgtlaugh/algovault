// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/zalgorithm
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/suffix_array.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;
    int n = s.size();
    SuffixArrayQueries<string> sa(s);

    for (int i = 0; i < n; i++) cout << sa.common_prefix(0, i) << (i + 1 < n ? ' ' : '\n');
    return 0;
}
