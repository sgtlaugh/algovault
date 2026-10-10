// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/lcm_convolution
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/gcd_lcm_convolution.cpp"
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
    int n = read_int();
    vector<long long> a(n + 1), b(n + 1);
    for (int i = 1; i <= n; i++) a[i] = read_int();
    for (int i = 1; i <= n; i++) b[i] = read_int();

    auto c = DivisorLattice(n).lcm_convolution(a, b, 998244353);
    print_line(vector<long long>(c.begin() + 1, c.end()));
    return 0;
}
