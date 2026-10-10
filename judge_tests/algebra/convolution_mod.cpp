// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/convolution_mod
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

long long read_int(){
    int c = getchar_unlocked();
    while (c != '-' && (c < '0' || c > '9')) c = getchar_unlocked();
    bool neg = c == '-';
    if (neg) c = getchar_unlocked();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return neg ? -x : x;
}

void write_int(long long x, char end){
    char s[24];
    int n = 0;
    if (x < 0) putchar_unlocked('-'), x = -x;
    do s[n++] = '0' + x % 10; while (x /= 10);
    while (n) putchar_unlocked(s[--n]);
    putchar_unlocked(end);
}

template <typename T>
void write_all(const vector<T>& v){
    for (size_t i = 0; i < v.size(); i++) write_int(v[i], i + 1 == v.size() ? '\n' : ' ');
}

int main(){
    int n = read_int(), m = read_int();
    vector<long long> a(n), b(m);
    for (auto& x : a) x = read_int();
    for (auto& x : b) x = read_int();

    write_all(ntt::multiply(a, b));
    return 0;
}
