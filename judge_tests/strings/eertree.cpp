// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/eertree
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/palindromic_tree.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    cin >> s;
    PalindromicTree<> tree;
    vector<int> longest_suffix;
    longest_suffix.reserve(s.size());
    for (char c : s){
        tree.add(c);
        longest_suffix.push_back(tree.last);
    }

    /// Node 0 (empty) is EVEN = 0, node 1 (imaginary) is ODD = -1, real nodes v >= 2 are created in judge order as v - 1
    auto id = [](int v){ return v < 2 ? -v : v - 1; };
    int nodes = tree.len.size();
    vector<int> parent(nodes);
    for (int u = 0; u < nodes; u++){
        for (int v : tree.next[u]){
            if (v != -1) parent[v] = u;
        }
    }

    cout << nodes - 2 << '\n';
    for (int v = 2; v < nodes; v++) cout << id(parent[v]) << ' ' << id(tree.link[v]) << '\n';
    for (int i = 0; i < (int)s.size(); i++) cout << id(longest_suffix[i]) << (i + 1 < (int)s.size() ? ' ' : '\n');
    return 0;
}
