// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/segment_tree_beats.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != '-' && (c < '0' || c > '9')) c = getchar();

    bool negative = c == '-';
    if (negative) c = getchar();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return negative ? -x : x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<long long> a(n);
    for (auto& x: a) x = read_int();

    SegmentTreeBeats st(a);
    while (q--){
        int t = read_int(), l = read_int() + 1, r = read_int();
        if (t == 3){
            printf("%lld\n", st.query_sum(l, r));
            continue;
        }

        long long b = read_int();
        if (t == 0) st.chmin(l, r, b);
        else if (t == 1) st.chmax(l, r, b);
        else st.add(l, r, b);
    }
    return 0;
}
