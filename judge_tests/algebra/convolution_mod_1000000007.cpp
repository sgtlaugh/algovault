// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod_1000000007
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c != '-' && !isdigit(c)) c = getchar();

    bool neg = c == '-';
    if (neg) c = getchar();

    long long x = 0;
    for (; isdigit(c); c = getchar()) x = x * 10 + (c - '0');
    return neg ? -x : x;
}

template<typename T>
void print_line(const vector<T>& v){
    string out;
    for (size_t i = 0; i < v.size(); i++){
        if (i) out += ' ';
        out += to_string(v[i]);
    }
    out += '\n';
    fwrite(out.data(), 1, out.size(), stdout);
}

int main(){
    int n = read_int(), m = read_int();
    vector<long long> a(n), b(m);
    for (auto& x : a) x = read_int();
    for (auto& x : b) x = read_int();

    print_line(ntt::mod_multiply(a, b, 1000000007));
    return 0;
}
