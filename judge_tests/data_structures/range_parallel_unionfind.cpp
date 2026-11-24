// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_parallel_unionfind
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/range_parallel_dsu.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    const long long MOD = 998244353;
    int n = read_int(), q = read_int();
    vector<long long> sum(n);
    for (auto& x : sum) x = read_int();

    /// Merging sets with weight sums s and t adds s * t to the sum over pairs in one set
    RangeParallelDSU dsu(n);
    long long total = 0;
    string out;
    while (q--){
        int k = read_int(), a = read_int(), b = read_int();
        dsu.connect(a, b, k, [&](int root, int absorbed){
            total = (total + sum[root] * sum[absorbed]) % MOD;
            sum[root] = (sum[root] + sum[absorbed]) % MOD;
        });
        out += to_string(total) + '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
