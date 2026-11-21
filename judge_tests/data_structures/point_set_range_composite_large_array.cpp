// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite_large_array
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/dynamic_segment_tree.cpp"
#undef main

const long long MOD = 998244353;

/// int fields keep a node at 28 bytes, so 5e5 point updates of about 30 nodes each fit the memory limit
struct Affine{
    int a = 1, b = 0;
};

/// combine(f, g) applies f first, then g, matching the left to right order of query
struct Composite{
    static Affine combine(Affine f, Affine g){ return {int(1LL * f.a * g.a % MOD), int((1LL * f.b * g.a + g.b) % MOD)}; }
    static Affine apply(Affine, Affine lz, long long){ return lz; }
    static Affine compose(Affine, Affine new_lz){ return new_lz; }
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

    DynamicSegmentTree<Affine, Composite> seg(0, n - 1);
    seg.nodes.reserve(32LL * q + 2);
    while (q--){
        int t = read_int(), x = read_int(), y = read_int(), z = read_int();
        if (t == 0) seg.update(x, x, Affine{y, z});
        else{
            Affine f = seg.query(x, y - 1);
            printf("%lld\n", (1LL * f.a * z + f.b) % MOD);
        }
    }

    return 0;
}
