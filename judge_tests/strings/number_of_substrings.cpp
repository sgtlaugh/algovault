// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/number_of_substrings
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/suffix_automaton.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;

    cout << SuffixAutomaton<>(s).distinct_substrings() << '\n';
    return 0;
}
