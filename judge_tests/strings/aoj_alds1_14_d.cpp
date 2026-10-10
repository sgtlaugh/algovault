// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_14_D
// competitive-verifier: TLE 3
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/suffix_array.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string t, p;
    int q;
    cin >> t >> q;
    SuffixArrayQueries<string> sa(t);

    while (q--){
        cin >> p;
        auto [lo, hi] = sa.occurrences(p);
        cout << (hi > lo) << '\n';
    }
    return 0;
}
