// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/DSL_2_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/kd_tree.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;

    vector<Point> points(n);
    for (auto& p : points){
        if (scanf("%lld %lld", &p.x, &p.y) != 2) return 0;
    }
    KDTree tree(points);

    int q;
    if (scanf("%d", &q) != 1) return 0;
    string out;
    while (q--){
        long long sx, tx, sy, ty;
        if (scanf("%lld %lld %lld %lld", &sx, &tx, &sy, &ty) != 4) return 0;

        vector<int> ids = tree.report(sx, sy, tx, ty);
        sort(ids.begin(), ids.end());
        for (int id : ids) out += to_string(id) + '\n';
        out += '\n';
    }

    fputs(out.c_str(), stdout);
    return 0;
}
