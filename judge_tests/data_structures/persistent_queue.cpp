// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/persistent_queue
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/persistent_segment_tree.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c != '-' && (c < '0' || c > '9')) c = getchar();

    bool negative = c == '-';
    if (negative) c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return negative ? -x : x;
}

int main(){
    int q = read_int();
    PersistentSegmentTree<long long, FirstMerge> arr(vector<long long>(q, 0));

    /// Queue S_i is (array version, head, tail) at index i + 1, the empty S_{-1} at index 0
    vector<int> version(q + 1, 0), head(q + 1, 0), tail(q + 1, 0);
    for (int i = 1; i <= q; i++){
        int type = read_int(), t = read_int() + 1;
        version[i] = version[t], head[i] = head[t], tail[i] = tail[t];

        if (type == 0) version[i] = arr.set(version[t], tail[i]++, read_int());
        else printf("%lld\n", arr.get(version[t], head[i]++));
    }
    return 0;
}
