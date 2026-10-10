// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_F
// competitive-verifier: ERROR 1e-5
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/circle.cpp"
#undef main

int main(){
    Point p;
    Circle c;
    if (scanf("%Lf %Lf %Lf %Lf %Lf", &p.x, &p.y, &c.center.x, &c.center.y, &c.r) != 5) return 0;

    vector<Point> pts = tangents_from_point(c, p);
    sort(pts.begin(), pts.end(), [](const Point& u, const Point& v){ return tie(u.x, u.y) < tie(v.x, v.y); });
    for (const Point& t : pts) printf("%.10Lf %.10Lf\n", t.x, t.y);
    return 0;
}
