// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_range_sum_large_array
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/dynamic_segment_tree.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    int a = 1, b = 0;
};

/// int sums and maps keep a node at 24 bytes, so about 4 log2(n) nodes per update fit the memory limit
struct AffineSum{
    static int combine(int x, int y){ return (x + y) % MOD; }
    static int apply(int val, Affine lz, long long len){ return (1LL * lz.a * val + 1LL * lz.b * (len % MOD)) % MOD; }
    static Affine compose(Affine old_lz, Affine new_lz){
        return {int(1LL * new_lz.a * old_lz.a % MOD), int((1LL * new_lz.a * old_lz.b + new_lz.b) % MOD)};
    }
};

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();

    DynamicSegmentTree<int, AffineSum, Affine> tree(0, n - 1);
    tree.nodes.reserve(4LL * 31 * q + 2);
    while (q--){
        int type = read_int(), l = read_int(), r = read_int() - 1;
        if (type == 0){
            int b = read_int(), c = read_int();
            tree.update(l, r, Affine{b, c});
        }
        else printf("%d\n", tree.query(l, r));
    }

    return 0;
}
