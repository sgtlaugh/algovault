// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/min_plus_convolution_convex_arbitrary
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/dp/monotone_minima.cpp"
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
    int n = read_int(), m = read_int();
    vector<long long> a(n), b(m);
    for (auto& x : a) x = read_int();
    for (auto& x : b) x = read_int();

    /// b[j] + a[k - j] is Monge for convex a; outside the valid j a penalty sloping away keeps it totally monotone
    const long long INF = 2e18;
    auto f = [&](int k, int j) -> long long{
        int lo = max(0, k - n + 1), hi = min(k, m - 1);
        if (j < lo) return INF + (lo - j);
        if (j > hi) return INF + (j - hi);
        return b[j] + a[k - j];
    };

    int rows = n + m - 1;
    vector<int> best = smawk(rows, m, f);

    string out;
    for (int k = 0; k < rows; k++){
        if (k) out += ' ';
        out += to_string(f(k, best[k]));
    }
    out += '\n';
    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
