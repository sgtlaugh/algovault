// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/stirling_number_of_the_first_kind
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/combinatorics/stirling_numbers.cpp"
#undef main

int main(){
    const long long mod = 998244353;
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<long long> c = stirling_first(n, mod);

    string out;
    char buf[16];
    for (int k = 0; k <= n; k++){
        long long s = (n - k) % 2 ? (mod - c[k]) % mod : c[k];
        out.append(buf, snprintf(buf, sizeof(buf), "%lld ", s));
    }
    out.back() = '\n';

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
