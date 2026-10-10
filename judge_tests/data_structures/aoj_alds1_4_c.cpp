// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/ALDS1_4_C
// competitive-verifier: TLE 2
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/trie.cpp"
#undef main

/// A, C, G, T become a..d so each node holds 4 children instead of the 20 that raw 'A'..'T' would need
string encode(const char* s){
    string res;
    for (; *s; s++) res += "abcd"[*s == 'A' ? 0 : *s == 'C' ? 1 : *s == 'G' ? 2 : 3];
    return res;
}

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    Trie<4, 'a'> trie;
    char command[8], word[16];
    while (n--){
        if (scanf("%7s %15s", command, word) != 2) return 0;
        if (command[0] == 'i') trie.insert(encode(word));
        else puts(trie.count_word(encode(word)) > 0 ? "yes" : "no");
    }
    return 0;
}
