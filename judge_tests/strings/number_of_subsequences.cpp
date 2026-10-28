// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/number_of_subsequences
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/subsequence_automaton.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    vector<int> a(n);
    for (int& x : a) cin >> x;

    cout << distinct_subsequences(a) << '\n';
    return 0;
}
