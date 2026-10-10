// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_14_C
// competitive-verifier: TLE 3
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/2D_pattern_matcher.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int h, w, r, c;
    cin >> h >> w;
    vector<string> text(h);
    for (auto& row : text) cin >> row;
    cin >> r >> c;
    vector<string> pattern(r);
    for (auto& row : pattern) cin >> row;

    auto matches = PatternMatcher2D(pattern).find(text);
    sort(matches.begin(), matches.end());

    for (auto [i, j] : matches) cout << i << ' ' << j << '\n';
    return 0;
}
