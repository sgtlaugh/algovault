// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_2_B
// competitive-verifier: TLE 1
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    while (q--){
        Point p[4];
        for (auto& pt : p){
            if (scanf("%lld %lld", &pt.x, &pt.y) != 2) return 0;
        }
        puts(segments_intersect(p[0], p[1], p[2], p[3]) ? "1" : "0");
    }
    return 0;
}
