// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/exp_of_formal_power_series
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/online_convolution.cpp"
#undef main

#define main ntt_main
#include "../../code_library/algebra/ntt.cpp"
#undef main

long long read_int(){
    int c = getchar_unlocked();
    while (c < '0' || c > '9') c = getchar_unlocked();

    long long x = 0;
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

void write_int(long long x, char end){
    char s[24];
    int n = 0;
    do s[n++] = '0' + x % 10; while (x /= 10);
    while (n) putchar_unlocked(s[--n]);
    putchar_unlocked(end);
}

int main(){
    const long long MOD = 998244353;
    int n = read_int();
    vector<long long> g(n), inv(n + 1, 1);
    for (int k = 0; k < n; k++) g[k] = read_int() * k % MOD;
    for (int i = 2; i <= n; i++) inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;

    /// f = e^a means f' = a' f, so i f_i = sum over j < i of f_j (i - j) a_(i-j)
    auto multiply = [](const vector<long long>& a, const vector<long long>& b){ return ntt::multiply(a, b); };
    auto next = [&](int i, long long s){ return i ? s * inv[i] % MOD : 1; };
    vector<long long> f = online_convolution(n, g, MOD, multiply, next);

    for (int i = 0; i < n; i++) write_int(f[i], i + 1 == n ? '\n' : ' ');
    return 0;
}
