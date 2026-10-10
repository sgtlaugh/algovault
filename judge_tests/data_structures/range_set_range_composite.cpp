// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_set_range_composite
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    long long a = 1, b = 0;
};

Affine then(Affine first, Affine second){
    return {second.a * first.a % MOD, (second.a * first.b + second.b) % MOD};
}

/// The library is customized by editing its CHANGE lines, a judge test specializes those members instead
template<>
Affine SegmentTree<Affine, Affine>::merge(Affine left, Affine right){
    return then(left, right);
}

/// A range of len positions all set to lz composes to lz applied len times
template<>
Affine SegmentTree<Affine, Affine>::apply(Affine, Affine lz, int len){
    Affine res;
    for (; len > 0; len >>= 1, lz = then(lz, lz)){
        if (len & 1) res = then(res, lz);
    }
    return res;
}

template<>
Affine SegmentTree<Affine, Affine>::compose(Affine, Affine new_lz){
    return new_lz;
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

    vector<Affine> f(n);
    for (Affine& g: f) g.a = read_int(), g.b = read_int();

    SegmentTree<Affine, Affine> tree(f, Affine());
    string out;
    while (q--){
        int type = read_int(), l = read_int() + 1, r = read_int();
        if (type == 0){
            long long c = read_int(), d = read_int();
            tree.update(l, r, {c, d});
        }
        else{
            long long x = read_int();
            Affine g = tree.query(l, r);
            out += to_string((g.a * x + g.b) % MOD) + '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
