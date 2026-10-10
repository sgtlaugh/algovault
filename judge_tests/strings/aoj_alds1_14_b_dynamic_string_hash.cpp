// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_14_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/dynamic_string_hash.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string t, p;
    cin >> t >> p;
    int n = t.size(), m = p.size();
    if (m > n) return 0;

    DynamicStringHash ht(t), hp(p);
    unsigned long long target = hp.hash(0, m - 1);

    for (int i = 0; i + m <= n; i++){
        if (ht.hash(i, i + m - 1) == target) cout << i << '\n';
    }
    return 0;
}
