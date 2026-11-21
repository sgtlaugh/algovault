// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/persistent_range_affine_range_sum
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/persistent_segment_tree.cpp"
#undef main

const long long MOD = 998244353;

struct Affine{
    long long a = 1, b = 0;
};

struct AffineSum{
    static long long combine(long long x, long long y){ return (x + y) % MOD; }
    static long long apply(long long val, Affine lz, int len){ return (lz.a * val + lz.b * len) % MOD; }
    static Affine compose(Affine old_lz, Affine new_lz){
        return {new_lz.a * old_lz.a % MOD, (new_lz.a * old_lz.b + new_lz.b) % MOD};
    }
};

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

    bool negative = c == '-';
    if (negative) c = read_char();
    int x = 0;
    for (; c >= '0' && c <= '9'; c = read_char()) x = x * 10 + (c - '0');
    return negative ? -x : x;
}

int main(){
    int n = read_int(), q = read_int();

    vector<long long> a(n);
    for (long long& x: a) x = read_int();

    PersistentLazySegmentTree<long long, AffineSum, Affine> seg(a);

    /// Query i's version is version_of[i + 1], A_{-1} is version 0; a sum query keeps the version it read
    vector<int> version_of(q + 1, 0);
    string out;
    for (int i = 0; i < q; i++){
        int type = read_int(), k = version_of[read_int() + 1];
        if (type == 0){
            int l = read_int(), r = read_int() - 1;
            long long b = read_int(), c = read_int();
            version_of[i + 1] = seg.update(k, l, r, {b, c});
        }
        else if (type == 1){
            int s = version_of[read_int() + 1], l = read_int(), r = read_int() - 1;
            version_of[i + 1] = seg.copy(k, s, l, r);
        }
        else{
            int l = read_int(), r = read_int() - 1;
            out += to_string(seg.query(k, l, r)) + '\n';
            version_of[i + 1] = k;
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
