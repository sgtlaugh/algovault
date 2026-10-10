// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_C
// competitive-verifier: TLE 1
// competitive-verifier: ERROR 1e-6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Point p[3];
    for (auto& q : p){
        long long x, y;
        if (scanf("%lld %lld", &x, &y) != 2) return 0;
        q = {(Float)x, (Float)y};
    }

    Circle c = circumcircle(p[0], p[1], p[2]);
    printf("%.12Lf %.12Lf %.12Lf\n", c.center.x, c.center.y, c.r);
    return 0;
}
