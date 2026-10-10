// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/jump_on_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/long_path_decomposition.cpp"
#undef main

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    LongPathDecomposition tree(n);
    for (int i = 0; i + 1 < n; i++){
        int a = read_int(), b = read_int();
        tree.add_edge(a, b);
    }
    tree.build(0);

    /// The library has no lca, so binary search the deepest depth d where both ancestors at depth d coincide
    auto lca = [&](int u, int v){
        int lo = 0, hi = min(tree.depth(u), tree.depth(v));
        while (lo < hi){
            int d = (lo + hi + 1) / 2;
            if (tree.kth_ancestor(u, tree.depth(u) - d) == tree.kth_ancestor(v, tree.depth(v) - d)) lo = d;
            else hi = d - 1;
        }
        return tree.kth_ancestor(u, tree.depth(u) - lo);
    };

    while (q--){
        int s = read_int(), t = read_int(), i = read_int();
        int l = lca(s, t), up_len = tree.depth(s) - tree.depth(l), len = up_len + tree.depth(t) - tree.depth(l);
        int res = i > len ? -1 : i <= up_len ? tree.kth_ancestor(s, i) : tree.kth_ancestor(t, len - i);
        printf("%d\n", res);
    }
    return 0;
}
