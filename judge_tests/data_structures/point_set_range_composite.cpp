// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/point_set_range_composite
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    long long a = 1, b = 0;
};

/// merge(f, g) applies f first, then g, matching the left to right order of query
template<>
Affine SegmentTree<Affine>::merge(Affine f, Affine g){
    return {f.a * g.a % MOD, (f.b * g.a + g.b) % MOD};
}

template<>
Affine SegmentTree<Affine>::apply(Affine, Affine lz, int){
    return lz;
}

template<>
Affine SegmentTree<Affine>::compose(Affine, Affine new_lz){
    return new_lz;
}

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<Affine> fs(n);
    for (auto& f: fs) f.a = read_int(), f.b = read_int();

    auto seg = SegmentTree<Affine>(fs, Affine{1, 0});
    while (q--){
        int t = read_int(), x = read_int(), y = read_int(), z = read_int();
        if (t == 0) seg.update(x + 1, x + 1, Affine{y, z});
        else{
            Affine f = seg.query(x + 1, y);
            printf("%lld\n", (f.a * z + f.b) % MOD);
        }
    }
    return 0;
}
