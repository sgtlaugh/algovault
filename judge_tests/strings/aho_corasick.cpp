// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/aho_corasick
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/strings/aho_corasick.cpp"
#undef main

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int q;
    cin >> q;
    AhoCorasick ac;
    vector<int> terminal(q);
    string s;
    for (int i = 0; i < q; i++){
        cin >> s;
        ac.insert(s);

        int cur = 0;
        for (char ch : s) cur = ac.go[cur][ac.edge[(unsigned char)ch]];
        terminal[i] = cur;
    }

    /// build() overwrites go with full transitions, so the trie parents must be read before it
    int n = ac.go.size();
    vector<int> parent(n, -1);
    for (int u = 0; u < n; u++){
        for (int v : ac.go[u]){
            if (v) parent[v] = u;
        }
    }
    ac.build();

    cout << n << '\n';
    for (int v = 1; v < n; v++) cout << parent[v] << ' ' << ac.fail[v] << '\n';
    for (int i = 0; i < q; i++) cout << terminal[i] << (i + 1 < q ? ' ' : '\n');
    return 0;
}
