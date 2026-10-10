// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_E
// competitive-verifier: TLE 1
// competitive-verifier: ERROR 1e-6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Circle c[2];
    for (auto& ci : c){
        long long x, y, r;
        if (scanf("%lld %lld %lld", &x, &y, &r) != 3) return 0;
        ci = {{(Float)x, (Float)y}, (Float)r};
    }

    vector<Point> pts = circle_circle_intersection(c[0], c[1]);
    if (pts.size() == 1) pts.push_back(pts[0]);

    /// Equal x values come out rounded apart, so a tie on x is decided by y within a tolerance
    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b){
        if (fabsl(a.x - b.x) > 1e-9L) return a.x < b.x;
        return a.y < b.y;
    });
    printf("%.12Lf %.12Lf %.12Lf %.12Lf\n", pts[0].x, pts[0].y, pts[1].x, pts[1].y);
    return 0;
}
