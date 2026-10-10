// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/kth_term_of_linearly_recurrent_sequence
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/algebra/linear_recurrence.cpp"
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

int main(){
    const int mod = 998244353;
    int d = read_int();
    long long k = read_int();

    vector<int> a(d), c(d);
    for (auto& x : a) x = read_int();
    for (auto& x : c) x = read_int();

    /// the constructor takes c in the order f(x) = sum c[i] * f(x - d + i), the reverse of c_1 .. c_d
    reverse(c.begin(), c.end());
    printf("%d\n", LinearRecurrence(a, mod, c).nth_term(k));
    return 0;
}
