// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_add_range_min
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree.cpp"
#undef main

/// A distinct value type, so specializing merge and apply leaves SegmentTree<long long> untouched
struct MinValue{
    long long v;
};

template<>
MinValue SegmentTree<MinValue, long long>::merge(MinValue a, MinValue b){
    return {min(a.v, b.v)};
}

template<>
MinValue SegmentTree<MinValue, long long>::apply(MinValue val, long long lz, int){
    return {val.v + lz};
}

static char buf[1 << 25];
int buf_len = 0, buf_pos = 0;

long long read_int(){
    while (buf_pos < buf_len && buf[buf_pos] != '-' && (buf[buf_pos] < '0' || buf[buf_pos] > '9')) buf_pos++;
    bool negative = buf_pos < buf_len && buf[buf_pos] == '-';
    if (negative) buf_pos++;

    long long x = 0;
    while (buf_pos < buf_len && buf[buf_pos] >= '0' && buf[buf_pos] <= '9') x = x * 10 + (buf[buf_pos++] - '0');
    return negative ? -x : x;
}

int main(){
    buf_len = fread(buf, 1, sizeof(buf), stdin);
    int n = read_int(), q = read_int();

    vector<MinValue> a(n);
    for (auto& x : a) x.v = read_int();
    auto seg = SegmentTree<MinValue, long long>(a, MinValue{LLONG_MAX});

    string out;
    while (q--){
        int t = read_int(), l = read_int(), r = read_int();
        if (t == 0) seg.update(l + 1, r, read_int());
        else{
            out += to_string(seg.query(l + 1, r).v);
            out += '\n';
        }
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
