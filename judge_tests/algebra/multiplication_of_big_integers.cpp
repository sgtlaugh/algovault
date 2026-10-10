// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/multiplication_of_big_integers
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

string read_token(){
    string s;
    int c = getchar();
    while (c != EOF && isspace(c)) c = getchar();
    while (c != EOF && !isspace(c)){
        s.push_back(c);
        c = getchar();
    }
    return s;
}

int main(){
    int t = stoi(read_token());

    string out;
    while (t--){
        Bignum a(read_token()), b(read_token());
        out += (a * b).to_string();
        out.push_back('\n');
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
