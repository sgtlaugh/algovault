// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/static_range_frequency
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/data_structures/coordinate_compression.cpp"
#undef main

int read_int(){
    int c = getchar();
    while (c < '0' || c > '9') c = getchar();

    int x = 0;
    for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
    return x;
}

int main(){
    int n = read_int(), q = read_int();
    vector<int> a(n);
    for (int& x : a) x = read_int();

    vector<int> vals = compress(a);
    vector<vector<int>> positions(vals.size());
    for (int i = 0; i < n; i++) positions[a[i]].push_back(i);

    string out;
    while (q--){
        int l = read_int(), r = read_int(), x = read_int();
        int id = lower_bound(vals.begin(), vals.end(), x) - vals.begin();
        int count = 0;
        if (id < (int)vals.size() && vals[id] == x){
            const vector<int>& p = positions[id];
            count = lower_bound(p.begin(), p.end(), r) - lower_bound(p.begin(), p.end(), l);
        }
        out += to_string(count) + '\n';
    }

    fwrite(out.data(), 1, out.size(), stdout);
    return 0;
}
