// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/2270
// competitive-verifier: TLE 3
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/tree_path_mo.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<int> x(n);
    for (int& v : x) v = read_int();

    TreePathMo tree(n);
    for (int i = 0; i < n - 1; i++){
        int a = read_int() - 1, b = read_int() - 1;
        tree.add_edge(a, b);
    }
    tree.build();

    vector<pair<int, int>> paths(q);
    vector<int> rank(q);
    for (int i = 0; i < q; i++) paths[i].first = read_int() - 1, paths[i].second = read_int() - 1, rank[i] = read_int();

    vector<int> vals = x;
    sort(vals.begin(), vals.end());
    vals.erase(unique(vals.begin(), vals.end()), vals.end());
    vector<int> id(n);
    for (int v = 0; v < n; v++) id[v] = lower_bound(vals.begin(), vals.end(), x[v]) - vals.begin() + 1;

    /// Counts of the window's values in a Fenwick tree; the l-th smallest is found by descending its powers of two
    int m = vals.size(), top = 1;
    while (top * 2 <= m) top *= 2;
    vector<int> fen(m + 1, 0);
    auto change = [&](int v, int d){ for (int i = id[v]; i <= m; i += i & -i) fen[i] += d; };
    auto kth = [&](int k){
        int pos = 0;
        for (int step = top; step; step >>= 1){
            if (pos + step <= m && fen[pos + step] < k) pos += step, k -= fen[pos];
        }

        return vals[pos];
    };

    vector<int> res = tree.run(paths, [&](int v){ change(v, 1); }, [&](int v){ change(v, -1); }, [&](int i){ return kth(rank[i]); });

    string out;
    for (int v : res) out += to_string(v) + '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
