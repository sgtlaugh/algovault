// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/rooted_tree_isomorphism_classification
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/tree_isomorphism.cpp"
#undef main

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int();
    IsoTree tree(n);
    for (int i = 1; i < n; i++) tree.add_edge(read_int(), i);

    /// A fresh canonizer hands out ids 0..K-1 in order, so the dictionary size is K
    TreeCanonizer canon;
    vector<int> id = canon.subtree_ids(tree, 0);

    string out = to_string(canon.ids.size()) + '\n';
    for (int v = 0; v < n; v++) out += to_string(id[v]) + (v + 1 < n ? ' ' : '\n');
    fputs(out.c_str(), stdout);
    return 0;
}
