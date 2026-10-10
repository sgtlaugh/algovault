// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_2_C
// competitive-verifier: TLE 1
// competitive-verifier: ERROR 1e-8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    int q;
    if (scanf("%d", &q) != 1) return 0;

    while (q--){
        PointF a, b, c, d;
        if (scanf("%lf %lf %lf %lf %lf %lf %lf %lf", &a.x, &a.y, &b.x, &b.y, &c.x, &c.y, &d.x, &d.y) != 8) return 0;
        PointF res = segment_intersection(a, b, c, d)[0];
        printf("%.10f %.10f\n", res.x, res.y);
    }
    return 0;
}
