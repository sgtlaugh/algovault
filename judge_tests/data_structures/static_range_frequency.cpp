// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_frequency
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/merge_sort_tree.cpp"
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
    vector<int> a(n);
    for (auto& x: a) x = read_int();

    WaveletMatrix<int> wavelet(a);
    while (q--){
        int l = read_int(), r = read_int(), x = read_int();
        if (l == r) puts("0");
        else printf("%d\n", wavelet.count_equal(l, r - 1, x));
    }
    return 0;
}
