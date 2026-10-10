// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/manhattanmst
// competitive-verifier: TLE 5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/graphs/manhattan_mst.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<pair<long long, long long>> points(n);
    for (auto& [x, y] : points){
        if (scanf("%lld %lld", &x, &y) != 2) return 0;
    }

    auto res = manhattan_mst(points);

    string out = to_string(res.weight) + "\n";
    for (auto [u, v] : res.edges) out += to_string(u) + " " + to_string(v) + "\n";
    fputs(out.c_str(), stdout);
    return 0;
}
