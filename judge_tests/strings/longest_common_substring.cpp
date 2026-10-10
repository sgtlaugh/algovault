// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/longest_common_substring
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/suffix_array.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s, t;
    cin >> s >> t;
    int n = s.size(), m = t.size();

    /// -1 is a unique separator, so no common prefix runs from S into T
    vector<int> joined(s.begin(), s.end());
    joined.push_back(-1);
    joined.insert(joined.end(), t.begin(), t.end());
    auto result = suffix_array(joined);

    int best = 0, from_s = 0, from_t = 0;
    for (int i = 0; i + 1 < n + m + 1; i++){
        int x = result.sa[i], y = result.sa[i + 1];
        if (x > y) swap(x, y);
        if (x < n && y > n && result.lcp[i] > best) best = result.lcp[i], from_s = x, from_t = y - n - 1;
    }

    if (best == 0) cout << "0 0 0 0\n";
    else cout << from_s << ' ' << from_s + best << ' ' << from_t << ' ' << from_t + best << '\n';
    return 0;
}
