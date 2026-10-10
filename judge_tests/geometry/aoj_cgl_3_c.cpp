// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_3_C
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int n, q;
    if (scanf("%d", &n) != 1) return 0;
    vector<Point> poly(n);
    for (auto& p : poly){
        if (scanf("%lld %lld", &p.x, &p.y) != 2) return 0;
    }

    if (scanf("%d", &q) != 1) return 0;
    while (q--){
        Point t;
        if (scanf("%lld %lld", &t.x, &t.y) != 2) return 0;
        printf("%d\n", point_in_polygon(poly, t) + 1);
    }
    return 0;
}
