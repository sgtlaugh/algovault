// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/frequency_table_of_tree_distance
// competitive-verifier: TLE 10
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/trees/centroid_decomposition.cpp"
#undef main
#define main fft_main
#include "../../code_library/algebra/fft.cpp"
#undef main

int read_int(){
    int c = getchar_unlocked(), x = 0;
    while (c < '0' || c > '9') c = getchar_unlocked();
    for (; c >= '0' && c <= '9'; c = getchar_unlocked()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int();
    CentroidDecomposition cd(n);
    for (int i = 0; i + 1 < n; i++){
        int u = read_int(), v = read_int();
        cd.add_edge(u, v);
    }

    /// Ordered pairs through c are the square of all depth counts minus each branch's own square
    vector<long long> freq(n, 0);
    cd.build([&](int, const CentroidDecomposition::Branches& branches){
        vector<long long> all = {1};
        for (auto& branch : branches){
            vector<long long> cnt;
            for (auto [v, d] : branch){
                if ((int)cnt.size() <= d) cnt.resize(d + 1, 0);
                if ((int)all.size() <= d) all.resize(d + 1, 0);
                cnt[d]++, all[d]++;
            }
            vector<long long> sq = fft::square(cnt);
            for (int d = 0; d < (int)sq.size(); d++) freq[d] -= sq[d];
        }
        vector<long long> sq = fft::square(all);
        for (int d = 0; d < (int)sq.size(); d++) freq[d] += sq[d];
    });

    string out;
    for (int d = 1; d < n; d++) out += to_string(freq[d] / 2) + (d + 1 < n ? ' ' : '\n');
    fputs(out.c_str(), stdout);
    return 0;
}
