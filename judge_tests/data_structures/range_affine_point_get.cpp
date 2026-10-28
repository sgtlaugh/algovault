// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_affine_point_get
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

int read_char(){
    if (buf_pos == buf_len){
        buf_len = fread(buf, 1, sizeof(buf), stdin);
        buf_pos = 0;
        if (buf_len <= 0) return -1;
    }

    return buf[buf_pos++];
}

long long read_int(){
    int c = read_char();
    while (c != '-' && (c < '0' || c > '9')) c = read_char();
    bool neg = c == '-';
    if (neg) c = read_char();
    long long x = 0;
    for (; c >= '0' && c <= '9'; c = read_char()) x = x * 10 + c - '0';
    return neg ? -x : x;
}

const long long MOD = 998244353;

struct Affine{
    long long b, c;
};

struct AffineSum{
    static long long combine(long long a, long long b){ return (a + b) % MOD; }
    static long long apply(long long val, Affine lz, int len){ return (lz.b * val + lz.c * len) % MOD; }
    static Affine compose(Affine old_lz, Affine new_lz){
        return {old_lz.b * new_lz.b % MOD, (old_lz.c * new_lz.b + new_lz.c) % MOD};
    }
};

int main(){
    int n = read_int(), q = read_int();
    vector<long long> ar(n);
    for (auto& x: ar) x = read_int();

    SegmentTree<long long, AffineSum, Affine> st(ar);
    string out;
    while (q--){
        int t = read_int();
        if (t == 0){
            int l = read_int(), r = read_int() - 1;
            Affine lz;
            lz.b = read_int(), lz.c = read_int();
            st.update(l, r, lz);
        }
        else{
            int i = read_int();
            out += to_string(st.query(i, i));
            out += '\n';
        }
    }
    fputs(out.c_str(), stdout);
    return 0;
}
