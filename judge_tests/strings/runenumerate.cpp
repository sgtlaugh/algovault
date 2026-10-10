// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/runenumerate
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/tandem_repeats.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;
    vector<array<int, 3>> triples;
    tandem_repeats(s, [&](int first, int last, int l){
        triples.push_back({l, first, last});
    });
    sort(triples.begin(), triples.end());

    /// A maximal block of square starts [a, b] with half length l is a maximal substring s[a, b + 2l) with period l,
    /// which is the run whose minimum period t divides l, so each run shows up once per multiple of t that fits
    vector<array<int, 3>> runs;
    for (int i = 0; i < (int)triples.size();){
        auto [l, a, b] = triples[i];
        for (i++; i < (int)triples.size() && triples[i][0] == l && triples[i][1] <= b + 1; i++) b = max(b, triples[i][2]);
        runs.push_back({a, b + 2 * l, l});
    }
    sort(runs.begin(), runs.end());

    vector<array<int, 3>> res;
    for (int i = 0; i < (int)runs.size(); i++){
        if (i == 0 || runs[i][0] != runs[i - 1][0] || runs[i][1] != runs[i - 1][1]) res.push_back({runs[i][2], runs[i][0], runs[i][1]});
    }
    sort(res.begin(), res.end());

    cout << res.size() << '\n';
    for (auto [t, l, r] : res) cout << t << ' ' << l << ' ' << r << '\n';
    return 0;
}
