// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_H
// competitive-verifier: TLE 1
// competitive-verifier: ERROR 1e-5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    int n;
    Float r;
    if (scanf("%d %Lf", &n, &r) != 2) return 0;

    vector<Point> poly(n);
    for (auto& p : poly){
        if (scanf("%Lf %Lf", &p.x, &p.y) != 2) return 0;
    }

    printf("%.10Lf\n", circle_polygon_area({{0, 0}, r}, poly));
    return 0;
}
