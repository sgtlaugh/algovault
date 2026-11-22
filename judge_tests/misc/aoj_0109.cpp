// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/0109
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/misc/expression_parser.cpp"
#undef main

int main(){
    int t;
    if (scanf("%d", &t) != 1) return 0;

    char buf[256];
    while (t--){
        if (scanf("%255s", buf) != 1) return 0;
        string s = buf;
        if (!s.empty() && s.back() == '=') s.pop_back();

        Expression e(s);
        printf("%lld\n", e.evaluate([](const string& token){ return stoll(token); }));
    }

    return 0;
}
