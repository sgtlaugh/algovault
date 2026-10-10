// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_D
// competitive-verifier: ERROR 1e-6
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Circle c;
    int q;
    if (scanf("%Lf %Lf %Lf %d", &c.center.x, &c.center.y, &c.r, &q) != 4) return 0;

    while (q--){
        Point a, b;
        if (scanf("%Lf %Lf %Lf %Lf", &a.x, &a.y, &b.x, &b.y) != 4) return 0;
        vector<Point> pts = circle_line_intersection(c, a, b);
        if (pts.size() == 1) pts.push_back(pts[0]);
        sort(pts.begin(), pts.end(), [](const Point& u, const Point& v){ return tie(u.x, u.y) < tie(v.x, v.y); });
        printf("%.10Lf %.10Lf %.10Lf %.10Lf\n", pts[0].x, pts[0].y, pts[1].x, pts[1].y);
    }
    return 0;
}
