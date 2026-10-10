// competitive-verifier: PROBLEM https://onlinejudge.u-aizu.ac.jp/problems/CGL_3_B
#include <bits/stdc++.h>

#define main library_main
#include "../../code_library/geometry/convex_hull.cpp"
#undef main

int main(){
    int n;
    if (scanf("%d", &n) != 1) return 0;
    vector<Point> polygon(n);
    for (auto& p : polygon){
        if (scanf("%" SCNd64 " %" SCNd64, &p.x, &p.y) != 2) return 0;
    }

    /// The judge data has an all-collinear "polygon" (outside its simple polygon constraint) and expects 1,
    /// while is_convex rejects degenerate input by contract
    bool collinear = true;
    for (const Point& p : polygon) collinear &= cross(polygon[0], polygon[1], p) == 0;
    puts(collinear || is_convex(polygon) ? "1" : "0");
    return 0;
}
