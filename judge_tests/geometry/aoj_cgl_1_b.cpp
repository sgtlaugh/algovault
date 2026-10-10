// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_1_B
// competitive-verifier: TLE 1
// competitive-verifier: ERROR 1e-8
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/geometry.cpp"
#undef main

int main(){
    PointF a, b;
    int q;
    if (scanf("%lf %lf %lf %lf %d", &a.x, &a.y, &b.x, &b.y, &q) != 5) return 0;

    while (q--){
        PointF p;
        if (scanf("%lf %lf", &p.x, &p.y) != 2) return 0;
        PointF r = reflect(p, a, b);
        printf("%.10f %.10f\n", r.x, r.y);
    }
    return 0;
}
