// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/set_xor_min
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/trie.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int q = read_int();
    BinaryTrie<30> trie;
    while (q--){
        int t = read_int(), x = read_int();
        if (t == 0 && trie.count(x) == 0) trie.insert(x);
        else if (t == 1 && trie.count(x) > 0) trie.erase(x);
        else if (t == 2) printf("%llu\n", trie.min_xor(x));
    }
    return 0;
}
