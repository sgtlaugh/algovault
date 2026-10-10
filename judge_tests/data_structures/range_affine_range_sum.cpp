// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_range_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

const long long MOD = 998244353;

struct Sum{
    long long value = 0;
};

struct Affine{
    long long a = 1, b = 0;
};

/// The library is customized by editing its CHANGE lines, a judge test specializes those members instead
template<>
Sum SegmentTree<Sum, Affine>::merge(Sum x, Sum y){
    return {(x.value + y.value) % MOD};
}

template<>
Sum SegmentTree<Sum, Affine>::apply(Sum val, Affine lz, int len){
    return {(lz.a * val.value + lz.b * len) % MOD};
}

template<>
Affine SegmentTree<Sum, Affine>::compose(Affine old_lz, Affine new_lz){
    return {new_lz.a * old_lz.a % MOD, (new_lz.a * old_lz.b + new_lz.b) % MOD};
}

static char buffer[1 << 16];
static size_t buffer_len = 0, buffer_pos = 0;

int read_char(){
    if (buffer_pos == buffer_len){
        buffer_len = fread(buffer, 1, sizeof(buffer), stdin);
        buffer_pos = 0;
        if (buffer_len == 0) return EOF;
    }
    return buffer[buffer_pos++];
}

int read_int(){
    int c = read_char();
    while (c == ' ' || c == '\n' || c == '\r') c = read_char();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = read_char()) x = x * 10 + (c - '0');
    return x;
}

int main(){
    int n = read_int(), q = read_int();

    vector<Sum> a(n);
    for (Sum& x: a) x.value = read_int();

    SegmentTree<Sum, Affine> tree(a);
    string out;
    while (q--){
        int type = read_int(), l = read_int() + 1, r = read_int();
        if (type == 0){
            long long b = read_int(), c = read_int();
            tree.update(l, r, {b, c});
        }
        else out += to_string(tree.query(l, r).value) + '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
