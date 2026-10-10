// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod_1000000007
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/number_theory/chinese_remainder_theorem.cpp"
#undef main
#define main ntt_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

long long read_int(){
    int c = getchar();
    while (c < '0') c = getchar();
    long long x = 0;
    for (; c >= '0'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

void write_int(long long x, char end){
    char buf[24];
    int len = 0;
    do buf[len++] = '0' + x % 10; while (x /= 10);
    while (len) putchar(buf[--len]);
    putchar(end);
}

int main(){
    int n = read_int(), m = read_int();
    vector<long long> a(n), b(m);
    for (auto& x : a) x = read_int();
    for (auto& x : b) x = read_int();

    auto c1 = ntt::multiply<ntt::P1, 3>(a, b), c2 = ntt::multiply<ntt::P2, 3>(a, b), c3 = ntt::multiply<ntt::P3, 11>(a, b);
    Garner g({ntt::P1, ntt::P2, ntt::P3}, 1000000007);

    vector<int64_t> rems(3);
    for (size_t i = 0; i < c1.size(); i++){
        rems[0] = c1[i], rems[1] = c2[i], rems[2] = c3[i];
        write_int(g(rems), i + 1 == c1.size() ? '\n' : ' ');
    }
    return 0;
}
