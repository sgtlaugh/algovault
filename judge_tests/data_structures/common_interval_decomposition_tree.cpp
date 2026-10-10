// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/common_interval_decomposition_tree
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/permutation_tree.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<int> p(n);
    for (int& x: p){
        if (scanf("%d", &x) != 1) return 0;
    }

    PermutationTree tree(p);

    /// The output needs parents numbered before children, so nodes are renumbered in preorder
    vector<int> id(tree.parent.size(), -1);
    vector<int> stack = {tree.root};
    string out;
    int count = 0;
    while (!stack.empty()){
        int v = stack.back();
        stack.pop_back();
        id[v] = count++;

        /// The checker labels leaves linear, as the reference solution does
        bool prime = tree.type[v] == PermutationTree::CUT && !tree.children[v].empty();
        int parent = tree.parent[v] == -1 ? -1 : id[tree.parent[v]];
        out += to_string(parent) + ' ' + to_string(tree.span[v].first) + ' ' + to_string(tree.span[v].second);
        out += prime ? " prime\n" : " linear\n";

        for (int child: tree.children[v]) stack.push_back(child);
    }

    printf("%d\n", count);
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
