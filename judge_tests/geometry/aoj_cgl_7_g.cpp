// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_7_G
// competitive-verifier: ERROR 1e-5
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

    vector<Point> pts;
    for (auto& [on_a, on_b] : common_tangents(c[0], c[1])) pts.push_back(on_a);

    /// Equal x values come out rounded apart, so a tie on x is decided by y within a tolerance
    sort(pts.begin(), pts.end(), [](const Point& a, const Point& b){
        if (fabsl(a.x - b.x) > 1e-9L) return a.x < b.x;
        return a.y < b.y;
    });
    for (const Point& p : pts) printf("%.12Lf %.12Lf\n", p.x, p.y);
    return 0;
}
