// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/enumerate_palindromes
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/hashing.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;
    int n = s.size();
    PolyHash h(vector<char>(s.begin(), s.end()));

    /// Center c spans s[c / 2 - k .. (c + 1) / 2 + k], binary search the largest k whose span reads the same backwards
    for (int c = 0; c < 2 * n - 1; c++){
        int left = c / 2, right = (c + 1) / 2, len = 0;
        if (s[left] == s[right]){
            int lo = 0, hi = min(left, n - 1 - right);
            while (lo < hi){
                int mid = (lo + hi + 1) / 2;
                if (h.get_hash(left - mid, right + mid) == h.rev_hash(left - mid, right + mid)) lo = mid;
                else hi = mid - 1;
            }
            len = 2 * lo + right - left + 1;
        }
        cout << len << (c + 1 < 2 * n - 1 ? ' ' : '\n');
    }
    return 0;
}
